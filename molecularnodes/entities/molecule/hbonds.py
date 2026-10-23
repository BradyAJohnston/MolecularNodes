import json
import logging
import math
import bpy
import databpy as db
import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
from ...blender import coll
from ...nodes import handlers as node_handlers
from ...nodes.geometry import StyleHBonds
from ..base import EntityType, MolecularEntity

logger = logging.getLogger(__name__)

# arguments to HydrogenBondAnalysis that are properties, adjustable in the UI
CUTOFFS = ("d_a_cutoff", "d_h_a_angle_cutoff", "d_h_cutoff")


class HBondCalculation:
    """
    Hydrogen bonds calculated with MDAnalysis, as the `hbonds` calculation

    Get it with `mol.calculations.hbonds()`, which sets it up if needed. Rather
    than an attribute, the calculation writes a data object (`object`) with one
    vertex per atom, in the same order as the atoms, and one edge per hydrogen
    bond between the hydrogen and the acceptor. The `Style HBonds` node draws
    the edges using positions sampled from the atoms. Uses
    `MDAnalysis.analysis.hydrogenbonds.HydrogenBondAnalysis`.
    """

    def __init__(self, entity: MolecularEntity):
        self._entity = entity
        self._hba = None
        self._hba_kwargs = None
        self._hba_trajectory = None
        self._kwargs = {}
        self._cutoffs = None
        self._setting_props = False
        self._uframe = None
        self._restore_failed = False

    @property
    def is_setup(self) -> bool:
        """Whether hydrogen bonds are set up as the `hbonds` calculation"""
        return self._hba is not None

    @property
    def object(self) -> bpy.types.Object | None:
        """The data object holding hydrogen bonds as edges, if created"""
        return bpy.data.objects.get(f".data_{self._entity.name}_hbonds")

    @property
    def _props(self):
        return self._entity.props.hbonds

    @property
    def cutoffs(self) -> dict[str, float]:
        """The cutoffs given to `HydrogenBondAnalysis`, from the properties

        `d_a_cutoff` and `d_h_cutoff` are in Angstroms, `d_h_a_angle_cutoff` in
        degrees. Change them with `mol.calculations.hbonds(d_a_cutoff=...)` or
        in the settings of `hbonds` in the Attributes panel.
        """
        props = self._props
        return {
            "d_a_cutoff": props.d_a_cutoff,
            "d_h_a_angle_cutoff": math.degrees(props.d_h_a_angle_cutoff),
            "d_h_cutoff": props.d_h_cutoff,
        }

    def _set_cutoffs(self, cutoffs: dict[str, float]) -> None:
        """Internal: Store cutoffs on the properties without recalculating"""
        values = {**cutoffs}
        values["d_h_a_angle_cutoff"] = math.radians(values["d_h_a_angle_cutoff"])
        self._setting_props = True
        try:
            for name, value in values.items():
                self._set_prop(name, value)
        finally:
            self._setting_props = False

    def _set_prop(self, name: str, value: any) -> None:
        """Internal: Set the property value only when different"""
        if getattr(self._props, name) != value:
            setattr(self._props, name, value)

    def _default_selections(self, selection: str) -> dict[str, str]:
        """Internal: Element based selections, as structure files have no charges to guess from"""
        universe = self._entity.universe
        # topologies such as GRO have no elements, so fall back to atom names.
        # Neutron structures model exchangeable hydrogens as deuterium
        if hasattr(universe.atoms, "elements"):
            h, n_o = "element H D", "element N O"
        else:
            h, n_o = "name H*", "name N* O*"
        # donors are paired with hydrogens by distance rather than bonds, as
        # hydrogens often have no bonds (deuterium names aren't in the CCD)
        return {
            "hydrogens_sel": f"({selection}) and ({h})",
            "acceptors_sel": f"({selection}) and ({n_o})",
            "donors_sel": f"({selection}) and ({n_o})",
        }

    def _setup_hba(self, selection: str, kwargs: dict) -> None:
        """Internal: Create the analysis for a selection and arguments"""
        universe = self._entity.universe
        kwargs = dict(kwargs)
        # cutoffs not given keep the values of the properties
        current = self.cutoffs
        cutoffs = {name: kwargs.pop(name, current[name]) for name in CUTOFFS}
        hba_kwargs = {
            "update_selections": False,
            **self._default_selections(selection),
            **cutoffs,
            **kwargs,
        }
        if not universe.select_atoms(hba_kwargs["hydrogens_sel"]):
            raise ValueError(
                "No hydrogens found, hydrogen bond analysis requires explicit hydrogens"
            )
        # create before assigning so a failing selection leaves the old state intact
        hba = HydrogenBondAnalysis(universe, **hba_kwargs)
        self._hba = hba
        self._hba_kwargs = hba_kwargs
        self._hba_trajectory = universe.trajectory
        self._kwargs = kwargs
        self._uframe = None
        self._set_prop("selection", selection)
        self._set_cutoffs(cutoffs)
        self._cutoffs = self.cutoffs
        try:
            self._set_prop("kwargs", json.dumps(kwargs))
        except TypeError:
            logger.warning("Hydrogen bond arguments can't be saved with the file")
            self._set_prop("kwargs", "{}")

    def _init(
        self,
        selection: str | None = None,
        updating: bool | None = None,
        **kwargs,
    ) -> None:
        """Internal: Set up hydrogen bonds if needed, and change the settings given

        See `CalculationManager.hbonds`, which calls this.
        """
        if self._entity._entity_type == EntityType.MD_STREAMING:
            raise ValueError("This does not apply for streamed trajectories")
        calculations = self._entity.calculations
        if self._hba is None:
            self._setup_hba("all" if selection is None else selection, kwargs)
            if self.object is None:
                db.create_object(
                    name=f".data_{self._entity.name}_hbonds",
                    collection=coll.data(),
                    vertices=np.zeros((len(self._entity.universe.atoms), 3)),
                )
            calculations.add(
                "hbonds",
                self._calculate,
                updating=True if updating is None else updating,
                refresh=False,
                data_object=self.object,
            )
        else:
            if (selection is not None and selection != self._props.selection) or kwargs:
                self._setup_hba(
                    self._props.selection if selection is None else selection,
                    {**self._kwargs, **kwargs},
                )
            if updating is not None:
                calculations.set_updating("hbonds", updating)
        # edges show straight away, whether or not they update on frame changes
        self._calculate(self._entity.universe)

    def _analysis(self) -> HydrogenBondAnalysis:
        """Internal: The analysis, rebuilt if the universe's trajectory reader
        was replaced (e.g. by `load_new` when the session is saved), as it reads
        frames through the reader it was created with"""
        universe = self._entity.universe
        if self._hba_trajectory is not universe.trajectory:
            self._hba = HydrogenBondAnalysis(universe, **self._hba_kwargs)
            self._hba_trajectory = universe.trajectory
        return self._hba

    def _update_from_props(self) -> None:
        """Internal: Recalculate after the cutoff properties change in the UI"""
        if self._hba is None or self._setting_props or self.cutoffs == self._cutoffs:
            return
        try:
            self._setup_hba(self._props.selection, self._kwargs)
        except Exception as e:
            logger.error(f"Hydrogen bonds failed with the new cutoffs: {e}")
            self._set_cutoffs(self._cutoffs)
            return
        self._calculate(self._entity.universe)

    def _restore(self) -> None:
        """Internal: Rebuild the `hbonds` calculation from its saved settings

        The calculation's function is not saved with a file, but its item on
        the object and the hydrogen bond properties are, so after a file is
        reopened it is set up again with the same selection, arguments and
        updating flag.
        """
        if self._hba is not None or self._restore_failed:
            return
        calculations = self._entity.calculations
        if "hbonds" in calculations:
            return
        item = calculations._item("hbonds")
        if item is None:
            return
        try:
            kwargs = json.loads(self._props.kwargs)
            self._init(self._props.selection, item.updating, **kwargs)
        except Exception as e:
            # don't retry on every frame change
            self._restore_failed = True
            logger.warning(f"Could not restore hydrogen bonds from saved settings: {e}")

    def _calculate(self, universe: mda.Universe) -> None:
        """Internal: Store hydrogen bonds of the current frame as edges"""
        trajectory = universe.trajectory
        uframe = trajectory.frame
        if self._hba is None or uframe == self._uframe:
            return None
        try:
            results = self._analysis().run(start=uframe, stop=uframe + 1).results
        finally:
            # running the analysis moves the trajectory
            if self._entity._entity_type != EntityType.MD_STREAMING:
                trajectory[uframe]
        # columns are frame, donor, hydrogen, acceptor, distance, angle
        pairs = results.hbonds[:, 2:4].astype(np.int32)
        # edges can't be removed from a mesh, so rebuild it each time
        mesh = self.object.data
        mesh.clear_geometry()
        mesh.vertices.add(len(universe.atoms))
        mesh.edges.add(len(pairs))
        mesh.edges.foreach_set("vertices", pairs.ravel())
        mesh.update()
        self._uframe = uframe
        return None


def _on_add_style_hbonds(
    node: bpy.types.GeometryNodeGroup, obj: bpy.types.Object | None
) -> None:
    """
    Setup for a ``Style HBonds`` node added from the asset library.

    Sets up hydrogen bonds for the molecule whose modifier tree received the
    node, if not already set up, and assigns the hydrogen bond data object to
    the node's Data Object input. A no-op when a data object is already
    assigned or the node isn't on a molecule.
    """
    from ...session import get_session
    from .base import Molecule

    socket = node.inputs.get("Data Object")
    if socket is None or socket.default_value is not None or obj is None:
        return
    entity = get_session().match(obj)
    if not isinstance(entity, Molecule):
        return
    try:
        hbonds = entity.calculations.hbonds()
    except ValueError as e:
        logger.warning(f"Hydrogen bonds not set up for {entity.name}: {e}")
        return
    socket.default_value = hbonds.object


node_handlers.register_on_add(StyleHBonds._name, _on_add_style_hbonds)
