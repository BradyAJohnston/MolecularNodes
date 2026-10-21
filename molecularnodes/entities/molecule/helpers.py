"""Helper classes for trajectory data management and Blender integration.

Provides utility classes for position caching, frame management, Blender property
synchronization, and attribute metadata collection.
"""

import logging
from typing import Optional
import databpy as db
import numpy as np
import numpy.typing as npt
from MDAnalysis import AtomGroup
from MDAnalysis.exceptions import NoDataError
from ...utils import (
    fraction,
    frame_mapper,
    frames_to_average,
    is_periodic_box,
    minimum_image_positions,
)

logger = logging.getLogger(__name__)


def _ag_to_bool(ag: AtomGroup) -> np.typing.NDArray[np.bool_]:
    """Convert AtomGroup to boolean mask for the entire universe."""
    mask = np.zeros(ag.universe.atoms.n_atoms, dtype=bool)
    mask[ag.ix] = True
    return mask


def _isin(values: np.ndarray, test_values) -> np.typing.NDArray[np.bool_]:
    """
    Per-element membership test, like `np.isin`, for arrays with few unique values.

    Residue and atom names are object arrays with many repeats, which `np.isin`
    compares one atom at a time against every test value. Testing only the unique
    values and mapping the result back is much faster for long lists of names.
    """
    unique, inverse = np.unique(values, return_inverse=True)
    test_values = set(test_values)
    is_in = np.fromiter((v in test_values for v in unique), bool, len(unique))
    return is_in[inverse].reshape(np.shape(values))


# ============================================================================
# Frame Management
# ============================================================================


# per-frame vector attributes read when the MDAnalysis Timestep has them:
# attribute name -> (Timestep attribute, whether to scale by world_scale).
# Velocities are lengths per time so they scale like positions; forces are not
# lengths and keep MDAnalysis units (kJ/(mol*Å))
TIMESTEP_ATTRIBUTES = {
    "velocity": ("velocities", True),
    "force": ("forces", False),
}

# names of the built-in attributes that per-atom `ts.data` arrays update: the
# B-factors and occupancies of each model of a PDB file, and LAMMPS dump
# charges. Other keys keep their name, or take a `ts_` prefix when it is the
# name of a built-in attribute
TS_DATA_NAMES = {"tempfactor": "b_factor", "occupancy": "occupancy", "q": "charge"}


class FrameManager:
    """Frame updates, position caching, interpolation, and periodic corrections.

    Encapsulates:
    - Caching of positions and other per-frame attributes
    - Frame-to-frame interpolation
    - Periodic boundary corrections
    - Frame averaging

    Maps between Blender scene frames and trajectory universe frames.

    Parameters
    ----------
    trajectory : Molecule
        Parent Molecule instance
    """

    def __init__(self, trajectory):
        """Initialize the frame manager.

        Parameters
        ----------
        trajectory : Molecule
            Parent Molecule instance
        """
        self.trajectory = trajectory
        # attribute name -> universe frame -> values, for positions and any
        # other per-frame attribute a frame has (see `_read_frame`)
        self.caches: dict[str, dict[int, np.ndarray]] = {}
        # universe frame -> unit cell dimensions (Å, degrees) of the frames in
        # the caches, or None for a frame without a box
        self._boxes: dict[int, Optional[np.ndarray]] = {}
        # attributes Molecular Nodes creates itself, which other `ts.data`
        # keys must not overwrite (`bond_type` is stored on the edges)
        self._builtin_names = {
            "position",
            "bond_type",
            *TIMESTEP_ATTRIBUTES,
            *getattr(trajectory, "_att_names", ()),
            *getattr(trajectory, "_blender_attributes", {}),
        }

    @property
    def attribute_caches(self) -> dict[str, dict[int, np.ndarray]]:
        """The caches of the per-frame attributes other than positions."""
        return {name: c for name, c in self.caches.items() if name != "position"}

    @property
    def n_frames(self) -> Optional[int]:
        """Get number of frames in trajectory.

        Returns None for streaming trajectories where the frame count
        is not known in advance.

        Returns
        -------
        Optional[int]
            Number of frames, or None for streaming trajectories
        """
        try:
            return self.trajectory.universe.trajectory.n_frames
        except (RuntimeError, AttributeError):
            # Streaming trajectories may not have known n_frames
            logger.debug("n_frames not available (likely streaming trajectory)")
            return None

    def _position_at_frame(self, frame: int) -> np.ndarray:
        """Get atom positions at a specific universe frame.

        Parameters
        ----------
        frame : int
            Universe frame number

        Returns
        -------
        np.ndarray
            Scaled atom positions
        """
        self.trajectory.uframe = frame
        return self.trajectory._scaled_position

    def _box_at_frame(self, frame: Optional[int]) -> Optional[np.ndarray]:
        """Unit cell dimensions of a universe frame, read when it was cached.

        Falls back to the box of the frame the universe is on for a frame that
        is not cached, or when `frame` is None.
        """
        if frame in self._boxes:
            return self._boxes[frame]
        return self.trajectory.universe.dimensions

    def adjust_periodic_positions(
        self, pos1: np.ndarray, pos2: np.ndarray, frame: Optional[int] = None
    ) -> np.ndarray:
        """Apply periodic boundary correction to positions.

        Shifts each atom of `pos2` to its periodic image nearest to it in
        `pos1` (the minimum image convention), so that interpolation and
        averaging take the short way across a boundary. Works for orthorhombic
        boxes and, via `MDAnalysis.lib.distances.minimize_vectors`, triclinic
        ones (see `minimum_image_positions`).

        Parameters
        ----------
        pos1 : np.ndarray
            Reference positions
        pos2 : np.ndarray
            Positions to correct
        frame : int, optional
            Universe frame whose box to use. When the box changes between
            frames (NPT), the box of the current frame is used.

        Returns
        -------
        np.ndarray
            Corrected positions (unchanged if correction is off or the frame
            has no usable box)
        """
        if not self.trajectory.correct_periodic:
            return pos2
        box = self._box_at_frame(frame)
        if not is_periodic_box(box):
            return pos2
        box = np.array(box, dtype=np.float64)
        box[:3] *= self.trajectory.world_scale
        return minimum_image_positions(pos1, pos2, box)

    def _frame_range(self, frame: int) -> npt.NDArray[np.int64]:
        """Get frame numbers to average over.

        Parameters
        ----------
        frame : int
            Center frame number

        Returns
        -------
        np.ndarray
            Frame numbers to include in average
        """
        n_frames = self.n_frames
        if n_frames is None:
            # Streaming trajectory - just return current frame
            return np.array([frame], dtype=np.int64)

        frames = frames_to_average(
            frame,
            n_frames,
            average=self.trajectory.average,
        )

        # Clamp frames to valid range
        frames = frames[(frames >= 0) & (frames < n_frames)]

        return frames

    def _frames_to_cache(
        self, frame: int, cache_ahead: bool = True
    ) -> npt.NDArray[np.int64]:
        """Get the frame numbers to hold in a cache for this frame.

        Extends a single-frame window by one look-ahead frame for
        interpolation, except at the end of the trajectory and for
        streaming trajectories (n_frames is None).

        Parameters
        ----------
        frame : int
            Current frame number
        cache_ahead : bool, default=True
            Whether to cache next frame for interpolation

        Returns
        -------
        np.ndarray
            Frame numbers to cache
        """
        frames_to_cache = self._frame_range(frame)
        n_frames = self.n_frames
        if (
            len(frames_to_cache) == 1
            and n_frames is not None
            and frames_to_cache[0] != (n_frames - 1)
            and cache_ahead
        ):
            frames_to_cache = np.array(
                (frames_to_cache[0], frames_to_cache[0] + 1), dtype=int
            )
        return frames_to_cache

    def _mapped_uframes(self, frame: int) -> tuple[int, int]:
        """Map a scene or manual frame to the current and next universe frames.

        Applies the subframe/offset mapping, clamps to trajectory bounds for
        non-streaming trajectories, and updates the UI frame property.

        The offset only applies when updating with the scene: it shifts where
        playback starts on the timeline, which has no meaning for a manually
        chosen frame. Subframes apply in both modes, so the manual frame is
        stepped through subframes the same way the scene frame is.

        Parameters
        ----------
        frame : int
            Scene frame number, or the manual frame when not updating with the
            scene

        Returns
        -------
        tuple[int, int]
            Current and next universe frame numbers
        """
        uframe_current = frame_mapper(
            frame=frame,
            subframes=self.trajectory.subframes,
            offset=self.trajectory.offset if self.trajectory.update_with_scene else 0,
        )
        uframe_next = uframe_current + 1

        n_frames = self.n_frames
        if n_frames is not None:
            last_frame = n_frames - 1
            if uframe_current >= last_frame:
                uframe_current = last_frame
                uframe_next = uframe_current

        # Update the frame_hidden property for the UI
        try:
            self.trajectory._frame = uframe_current
        except AttributeError:
            # Silently ignore if we can't write to Blender property
            pass

        return uframe_current, uframe_next

    def _trim_cache(self, frames_to_keep) -> None:
        """Drop cached frames other than these, once more than 2 are held."""
        if len(self.caches.get("position", {})) > 2:
            for cache in (*self.caches.values(), self._boxes):
                for f in [f for f in cache if f not in frames_to_keep]:
                    del cache[f]

    def _read_frame(self, frame: int) -> dict[str, np.ndarray]:
        """Read every per-frame attribute that a universe frame has.

        Positions are always read. The entity's `_att_names` (e.g. oxDNA base
        vectors) are read when the frame's `ts.data` has them, and velocities
        and forces when the Timestep has them: most trajectories have neither,
        and a TRR file can have them on only some frames. Other per-atom arrays
        in `ts.data` are read under their MN names (see `TS_DATA_NAMES`).

        Parameters
        ----------
        frame : int
            Universe frame number

        Returns
        -------
        dict[str, np.ndarray]
            Values by attribute name, in Blender units where they are lengths
        """
        values = {"position": self._position_at_frame(frame)}
        ts = self.trajectory.universe.trajectory.ts
        dimensions = ts.dimensions
        self._boxes[frame] = None if dimensions is None else np.array(dimensions)
        world_scale = self.trajectory.world_scale
        att_names = getattr(self.trajectory, "_att_names", ())
        for name in att_names:
            if name in ts.data:
                values[name] = ts.data[name] * world_scale
        for name, (ts_attribute, scaled) in TIMESTEP_ATTRIBUTES.items():
            if name in values or not getattr(ts, f"has_{ts_attribute}", False):
                continue
            try:
                data = getattr(ts, ts_attribute)
            except NoDataError:
                continue
            values[name] = data * world_scale if scaled else data.copy()
        # any other per-atom array the reader put in `ts.data`, such as PDB
        # occupancies and B-factors per model, LAMMPS dump `additional_columns`
        # or TNG `additional_blocks`. Single-frame structures keep the values
        # from import
        if self.n_frames != 1:
            for key, data in ts.data.items():
                name = TS_DATA_NAMES.get(key)
                if name is None:
                    name = f"ts_{key}" if key in self._builtin_names else key
                if key in att_names or name in values:
                    continue
                if isinstance(data, np.ndarray) and data.shape[:1] == (ts.n_atoms,):
                    values[name] = data.copy()
        return values

    def read_and_cache(self, frame: int) -> dict[str, np.ndarray]:
        """Read a universe frame afresh and cache it.

        Parameters
        ----------
        frame : int
            Universe frame number

        Returns
        -------
        dict[str, np.ndarray]
            Values by attribute name
        """
        values = self._read_frame(frame)
        self._trim_cache([frame])
        for name, value in values.items():
            self.caches.setdefault(name, {})[frame] = value
        return values

    def update_cache(self, frame: int, cache_ahead: bool = True) -> None:
        """Update the cache for the current frame.

        Intelligently caches based on averaging and interpolation settings,
        prefetching the next frame when needed.

        Parameters
        ----------
        frame : int
            Current frame number
        cache_ahead : bool, default=True
            Whether to cache next frame for interpolation
        """
        frames_to_cache = self._frames_to_cache(frame, cache_ahead)
        self._trim_cache(frames_to_cache)

        # Update the cache with any frames that are not yet cached
        for f in frames_to_cache:
            if f not in self.caches.get("position", {}):
                for name, value in self._read_frame(f).items():
                    self.caches.setdefault(name, {})[f] = value

    def cache_mean(self, frame: int) -> dict[str, np.ndarray]:
        """Get mean values from cached frames.

        Computes the average over the averaging window, with periodic
        corrections for positions if needed. Each attribute is averaged over
        the frames of the window that have it.

        Parameters
        ----------
        frame : int
            Center frame number

        Returns
        -------
        dict[str, np.ndarray]
            Mean values (or single frame if not averaging) of the attributes
            this frame has, by name
        """
        self.update_cache(frame)
        current = {
            name: cache[frame] for name, cache in self.caches.items() if frame in cache
        }
        if self.trajectory.average == 0:
            return current

        # park the universe on this frame so the box dimensions used by the
        # corrections below are a deterministic function of `frame` — after
        # cache updates the universe sits on whichever frame happened to be a
        # cache miss, and in an NPT box a borderline (half-box) correction
        # would flip between evaluations of the same frame
        self.trajectory.uframe = frame

        # average exactly the window for this frame — the cache can hold other
        # frames (interpolation look-ahead, a previous window), which must not
        # leak into the mean
        window = self._frame_range(frame)
        means = {}
        for name in current:
            if current[name].dtype.kind != "f":
                # integer and boolean values are not averaged
                means[name] = current[name]
                continue
            cache = self.caches[name]
            array = np.stack([cache[f] for f in window if f in cache])
            if name == "position" and self.trajectory.correct_periodic:
                # correct each frame against its (corrected) predecessor:
                # consecutive frames are where the nearest-image assumption
                # holds. Correcting everything against the window's first frame
                # mistakes real motion greater than half a box — legitimate in
                # unwrapped or per-molecule-imaged trajectories — for a boundary
                # crossing, and the windows of neighbouring frames then disagree
                # by a whole box
                for i in range(1, len(array)):
                    array[i] = self.adjust_periodic_positions(
                        array[i - 1], array[i], frame
                    )
            means[name] = np.mean(array, axis=0)
        return means

    def position_cache_mean(self, frame: int) -> np.ndarray:
        """Get mean positions from cached frames, see `cache_mean`."""
        return self.cache_mean(frame)["position"]

    def get_values_at_frame(self, frame: int) -> dict[str, np.ndarray]:
        """Get per-frame attributes for a given frame with all processing applied.

        Main entry point for frame updates. Handles frame mapping, interpolation,
        averaging, and periodic boundary corrections (for positions only). An
        attribute the current frame lacks is left out; one the next frame lacks
        is not interpolated.

        Parameters
        ----------
        frame : int
            Scene frame number, or the manual frame when not updating with the
            scene

        Returns
        -------
        dict[str, np.ndarray]
            Processed values by attribute name, always including `position`
        """
        uframe_current, uframe_next = self._mapped_uframes(frame)

        if self.trajectory.subframes > 0 and self.trajectory.interpolate:
            # Interpolate between current and next frame
            current = self.cache_mean(uframe_current)
            after = self.cache_mean(uframe_next)
            t = fraction(frame, self.trajectory.subframes + 1)
            values = {}
            for name, value in current.items():
                # integer and boolean values are not interpolated
                if name not in after or value.dtype.kind != "f":
                    values[name] = value
                    continue
                value_next = after[name]
                if name == "position":
                    # interpolate towards the nearest periodic image. Averaging
                    # only corrects within each window, so two averaged
                    # endpoints can still sit on opposite sides of the box — the
                    # correction applies to them just the same
                    value_next = self.adjust_periodic_positions(
                        value, value_next, uframe_current
                    )
                values[name] = db.lerp(value, value_next, t=t)
            return values
        elif self.trajectory.average > 0:
            # Return mean values for cached frames
            return self.cache_mean(uframe_current)
        else:
            # Just return current values
            return self.read_and_cache(uframe_current)

    def get_positions_at_frame(self, frame: int) -> np.ndarray:
        """Get positions for a given frame, see `get_values_at_frame`."""
        return self.get_values_at_frame(frame)["position"]
