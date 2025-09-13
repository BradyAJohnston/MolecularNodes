from __future__ import annotations
import json
import re
import urllib.request
from pathlib import Path
from typing import Any
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


def _json_get(url: str) -> Any:
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        charset = resp.headers.get_content_charset() or "utf-8"
        data = resp.read().decode(charset)
    return json.loads(data)


def _find_first_map_url(obj: Any, emdb_id: str) -> str | None:
    """Recursively search JSON for a .map.gz URL, prefer one containing the id."""
    found: list[str] = []
    if isinstance(obj, dict):
        for v in obj.values():
            url = _find_first_map_url(v, emdb_id)
            if url:
                found.append(url)
    elif isinstance(obj, list):
        for v in obj:
            url = _find_first_map_url(v, emdb_id)
            if url:
                found.append(url)
    elif isinstance(obj, str):
        if obj.endswith(".map.gz") and obj.startswith("http"):
            found.append(obj)

    if not found:
        return None
    # Prefer URLs that include the emdb id and standard emd_#### name
    for url in found:
        if f"emd_{emdb_id}" in url.lower() or f"EMD-{emdb_id}" in url:
            return url
    return found[0]


def fetch_emdb_map(emdb_id: str | int, progress: bool = False) -> Path:
    """Download (and cache) the main .map.gz for an EMDB entry using pooch.

    Tries the EMDB REST API to discover the map URL; if discovery fails,
    falls back to the well-known FTP URL pattern.
    """
    eid = _normalize_emdb_id(emdb_id)

    # Try to discover the file URL via REST API
    map_url: str | None = None
    api_candidates = [
        f"https://www.ebi.ac.uk/emdb/api/entry/EMD-{eid}",
        f"https://www.ebi.ac.uk/emdb/api/entry/EMD-{eid}/files",
    ]
    for api_url in api_candidates:
        try:
            payload = _json_get(api_url)
            candidate = _find_first_map_url(payload, eid)
            if candidate:
                map_url = candidate
                break
        except Exception:
            # Ignore and try next candidate or fallback
            continue

    # Fallback: direct FTP path pattern
    if not map_url:
        map_url = f"https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-{eid}/map/emd_{eid}.map.gz"

    cache_dir = pooch.os_cache("molecularnodes-tests")
    fname = f"emd_{eid}.map.gz"
    path = pooch.retrieve(
        url=map_url,
        known_hash=None,  # No checksum available a priori
        fname=fname,
        path=cache_dir,
        progressbar=progress,
    )
    return Path(path)
