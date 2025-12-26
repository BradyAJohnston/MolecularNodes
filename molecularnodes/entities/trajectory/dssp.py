import logging
import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis.dssp import DSSP, translate
from ..base import EntityType, MolecularEntity

logger = logging.getLogger(__name__)


class DSSPManager:
    """
    DSSP Manager for trajectories

    Show secondary structures computed using MDAnalysis.analysis.dssp
    """

    def __init__(self, entity: MolecularEntity):
        self._entity = entity
        self._DSSP = None
        self._dssp_results = None
        self._trajectory_average = None
        self._dssp_resindices = None
        self._dssp_vmap = np.vectorize({"H": 1, "E": 2, "-": 3}.get)
        self._props = None
        self._display_option = "none"
        self._window_size = 5
        self._threshold = None
        self._no_sec_struct = None

    def _set_dssp_resindices(self, resids: list) -> None:
        """Internal: Set resindices for DSSP resids"""
        if self._dssp_resindices is None:
            universe = self._entity.universe
            mask = np.isin(universe.residues.resids, resids)
            self._dssp_resindices = universe.residues.resindices[mask]

    def _get_sliding_window_indices(self, index: int) -> None:
        """Internal: Get sliding window indices based on current index"""
        half = self._window_size // 2
        indices = np.arange(index - half, index + half + 1)
        n = self._entity.universe.trajectory.n_frames
        mask = (indices >= 0) & (indices < n)
        return indices[mask]

    def _run_dssp(self, window: list) -> DSSP:
        """Internal: Run DSSP for a subset of frames"""
        trajectory = self._entity.universe.trajectory
        # save current frame
        current_frame = trajectory.frame
        try:
            dssp_results = self._DSSP.run(frames=window).results
        except Exception as e:
            dssp_results = None
            logger.debug(f"Failed to run DSSP for frames {window}: {e}")
        # restore current frame
        if self._entity._entity_type != EntityType.MD_STREAMING:
            trajectory[current_frame]
        return dssp_results

    def _calculate_sec_struct(self, universe: mda.Universe) -> np.ndarray:
        """Internal: Calculate secondary structures on frame change event"""
        if self._display_option == "none" or self._DSSP is None:
            return self._no_sec_struct
        # run dssp if required
        dssp_results = self._dssp_results
        frame = universe.trajectory.frame
        if self._display_option == "per-frame":
            if self._trajectory_average is None:
                dssp_results = self._run_dssp([frame])
        elif self._display_option == "sliding-window-average":
            frames = self._get_sliding_window_indices(frame)
            dssp_results = self._run_dssp(frames)
        # ensure we have dssp run results (per-frame or sliding window or full)
        if dssp_results is None:
            return self._no_sec_struct
        # set resindices corresponding to dssp resids
        self._set_dssp_resindices(dssp_results.resids)
        # create attribute data
        if self._display_option == "trajectory-average":
            dssp_chars = self._trajectory_average
        elif self._display_option == "sliding-window-average":
            dssp_chars = translate(dssp_results.dssp_ndarray.mean(axis=0))
        else:  # per-frame
            index = 0 if self._trajectory_average is None else frame
            dssp_chars = dssp_results.dssp[index]
        dssp_ints = self._dssp_vmap(dssp_chars)
        attribute_data = np.zeros(len(universe.atoms), dtype=int)
        attribute_data[self._dssp_resindices] = dssp_ints
        return attribute_data[universe.atoms.resindices]

    @property
    def _display_option_prop(self):
        """Internal: Getter for display option property"""
        if self._entity._entity_type == EntityType.MD_STREAMING:
            return self._props.display_option_streaming
        else:
            return self._props.display_option

    @_display_option_prop.setter
    def _display_option_prop(self, value: str) -> None:
        """Internal: Setter for display option property"""
        if self._entity._entity_type == EntityType.MD_STREAMING:
            self._props.display_option_streaming = value
        else:
            self._props.display_option = value

    def _ensure_init(self) -> None:
        """Internal: Ensure DSSP is initialized"""
        if self._DSSP is None:
            raise ValueError("DSSP is not initialized")

    def _ensure_no_streaming(self) -> None:
        """Internal: Ensure non streaming trajectory"""
        if self._entity._entity_type == EntityType.MD_STREAMING:
            raise ValueError("This does not apply for streamed trajectories")

    def _set_display_option(self, option: str) -> None:
        """Internal: Set the display option"""
        self._display_option = option
        if self._display_option_prop != option:
            self._display_option_prop = option

    def init(self) -> None:
        """
        Initialize DSSP
        """
        if self._DSSP is not None:
            raise ValueError("DSSP already initialized")
        universe = self._entity.universe
        self._DSSP = DSSP(universe)
        self._entity.calculations["sec_struct"] = self._calculate_sec_struct
        self._props = self._entity.object.mn.dssp
        # calculate no secondary structs attribute
        # protein - 3 (loop), rest - 0 (none)
        protein_resindices = universe.select_atoms("protein").residues.resindices
        no_sec_struct = np.zeros(len(universe.atoms), dtype=int)
        no_sec_struct[protein_resindices] = 3
        self._no_sec_struct = no_sec_struct[universe.atoms.resindices]
        # set and apply default
        self._set_display_option("none")
        self._props.applied = True

    def show_none(self) -> None:
        """
        Do not show secondary structures
        """
        self._ensure_init()
        self._set_display_option("none")

    def show_per_frame(self) -> None:
        """
        Show secondary structures calculated per frame
        """
        self._ensure_init()
        self._set_display_option("per-frame")

    def show_sliding_window_average(self, window_size: int = 5) -> None:
        """
        Show average secondary structures of a sliding window of frames

        Parameters
        ----------
        window_size: int, optional
            Size of the sliding window, default is 5 frames
        """
        self._ensure_init()
        self._ensure_no_streaming()
        self._window_size = window_size
        self._set_display_option("sliding-window-average")
        self._props.window_size = window_size
        self._props.applied = True

    def show_trajectory_average(self, threshold: float | None = None) -> None:
        """
        Show average secondary structures across all frames

        Parameters
        ----------
        threshold: float, optional
            Threshold to compare the mean against [0.0 - 1.0].
            When None, no threshold comparison is made
        """
        self._ensure_init()
        self._ensure_no_streaming()
        self._threshold = threshold
        if self._dssp_results is None:
            self._dssp_results = self._DSSP.run().results.copy()
        if threshold is not None:
            self._trajectory_average = translate(
                self._dssp_results.dssp_ndarray.mean(axis=0) > threshold
            )
            self._props.threshold = threshold
        else:
            self._trajectory_average = translate(
                self._dssp_results.dssp_ndarray.mean(axis=0)
            )
        self._set_dssp_resindices(self._dssp_results.resids)
        self._set_display_option("trajectory-average")
        self._props.applied = True
