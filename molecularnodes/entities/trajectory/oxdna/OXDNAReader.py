import numpy as np

from MDAnalysis.coordinates.base import ReaderBase
from MDAnalysis.lib import util


def _is_info_line(line: str):
    """
    Check if a line contains simulation information.

    Parameters
    ----------
    line : str
        Line from the OXDNA trajectory file.

    Returns
    -------
    bool
        True if line contains simulation information, False otherwise.
    """
    return line.startswith("t = ") or line.startswith("b = ") or line.startswith("E = ")


class OXDNAReader(ReaderBase):
    """
    Reader for OXDNA trajectory files.

    Parameters
    ----------
    filename : str
        Path to the OXDNA trajectory file.
    **kwargs : dict
        Additional arguments to pass to the reader, must include 'n_atoms'.

    Attributes
    ----------
    n_atoms : int
        Number of atoms in the system.
    ts : Timestep
        Current timestep object.
    n_frames : int
        Total number of frames in trajectory.
    """

    def __init__(self, filename, **kwargs):
        super(OXDNAReader, self).__init__(filename, **kwargs)

        self.n_atoms = kwargs["n_atoms"]
        self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)

        oxdnafile = self._oxdnafile = util.anyopen(filename, "rb")

        starts = []
        stops = []

        line = "a new line"
        previous_is_info = False
        current_is_info = False
        while line:
            line = oxdnafile.readline().decode()
            current_is_info = _is_info_line(line)

            # get the lines where the data values start and stop, which
            # are separated by the several lines of overall information
            if not current_is_info and previous_is_info:
                starts.append(oxdnafile.tell() - len(line))
            if current_is_info and not previous_is_info:
                stops.append(oxdnafile.tell() - len(line))

            previous_is_info = current_is_info

        # add the end of the file
        stops += [oxdnafile.tell()]

        self._start_offsets = starts
        self._stop_offsets = stops[1:]  # drop the first
        self.n_frames = len(self._start_offsets)

        self._read_frame(0)

    def _reopen(self):
        """
        Reopen the trajectory file and reset the timestep frame.
        """
        self.close()
        self._oxdnafile = util.anyopen(self.filename, "rb")
        self.ts.frame = -1

    def _read_next_timestep(self, ts=None):
        """
        Read the next timestep from the trajectory.

        Parameters
        ----------
        ts : Timestep, optional
            If provided, update this timestep object instead of creating a new one.

        Returns
        -------
        Timestep
            The timestep object containing the frame data.
        """
        frame = self.frame + 1
        return self._read_frame(frame)

    def _read_frame(self, frame):
        """
        Read a specific frame from the trajectory.

        Parameters
        ----------
        frame : int
            Frame number to read.

        Returns
        -------
        Timestep
            The timestep object containing the frame data.

        Raises
        ------
        OSError
            If the frame number is out of range.
        """
        try:
            start = self._start_offsets[frame]
            stop = self._stop_offsets[frame]
        except IndexError:
            raise OSError from None

        self._oxdnafile.seek(start)
        chunk = self._oxdnafile.read(stop - start)

        array = np.array(
            [
                np.array(line.split(), dtype=float)
                for line in chunk.decode().splitlines()
            ]
        )

        # TODO: also access and update the other values
        self.ts.positions = array[:, :3]

        for i, name in enumerate(
            ("base_vector", "base_normal", "velocity", "angular_velocity")
        ):
            starting_column = 3 * (i + 1)
            if starting_column >= array.shape[1]:
                continue
            cols = np.arange(3) + starting_column
            self.ts.data[name] = array[:, cols]
        self.ts.frame = frame

        return self.ts

    def close(self):
        """
        Close the trajectory file.
        """
        self._oxdnafile.close()
