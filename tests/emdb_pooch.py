from __future__ import annotations
import re
from pathlib import Path
import pooch


def _normalize_emdb_id(emdb_id: str | int) -> str:
    s = str(emdb_id).strip()
    s = s.upper()
    s = s.replace("EMD-", "").replace("EMD", "")
    # keep only digits
    s = re.sub(r"\D+", "", s)
    if not s:
        raise ValueError(f"Invalid EMDB id: {emdb_id}")
    return s


def fetch_emdb_map(emdb_id: str | int, progress: bool = False) -> Path:
    """Download (and cache) the main .map.gz for an EMDB entry using pooch.

    The archive is mirrored by every wwPDB site, so each is tried in turn: EBI's
    own host is often slow enough to time out from CI runners.
    """
    eid = _normalize_emdb_id(emdb_id)
    path = f"pub/emdb/structures/EMD-{eid}/map/emd_{eid}.map.gz"
    urls = [
        f"https://files.wwpdb.org/{path}",
        f"https://files.pdbj.org/{path}",
        f"https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-{eid}/map/emd_{eid}.map.gz",
    ]

    # give up on a host that doesn't answer quickly rather than waiting out the
    # CI step's timeout on it
    downloader = pooch.HTTPDownloader(progressbar=progress, timeout=(10, 60))
    errors = []
    for url in urls:
        try:
            return Path(
                pooch.retrieve(
                    url=url,
                    known_hash=None,  # No checksum available a priori
                    fname=f"emd_{eid}.map.gz",
                    path=pooch.os_cache("molecularnodes-tests"),
                    # pooch's Downloader protocol names the first argument
                    # differently from its own HTTPDownloader
                    downloader=downloader,  # ty: ignore[invalid-argument-type]
                )
            )
        except Exception as e:
            errors.append(f"{url}: {e}")
    raise ConnectionError("No EMDB mirror answered:\n" + "\n".join(errors))
