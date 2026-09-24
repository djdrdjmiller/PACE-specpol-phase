"""Explicit local path configuration for analysis notebooks."""

import os
from pathlib import Path
import tomllib


def load_paths(config_path=None):
    """Load paths from an explicit file, PACE_CONFIG, or an ancestor repo config.

    Values are expanded and resolved relative to the configuration file. No
    directories are created and no data are moved. Select a distinct config for
    each experiment; existing cache settings remain authoritative.
    """
    selected = config_path or os.environ.get("PACE_CONFIG")
    if selected:
        config = Path(selected).expanduser().resolve()
    else:
        config = next(
            (
                p / "config" / "paths.local.toml"
                for p in (Path.cwd(), *Path.cwd().parents)
                if (p / "config" / "paths.local.toml").is_file()
            ),
            None,
        )
        if config is None:
            raise FileNotFoundError(
                "Create config/paths.local.toml from paths.example.toml, or set "
                "PACE_CONFIG to an absolute local TOML configuration path. "
                "Use arcsix.example.toml for the Arctic notebook."
            )
    with config.open("rb") as stream:
        settings = tomllib.load(stream).get("paths", {})
    required = {"lut_root", "pair_cache", "work_root", "airborne_root", "export_root"}
    missing = required - settings.keys()
    if missing:
        raise ValueError(f"Missing paths in {config}: {', '.join(sorted(missing))}")
    result = {}
    for key in required:
        value = settings[key]
        if not isinstance(value, str) or not value.strip():
            raise ValueError(f"Path {key} must be a nonempty string in {config}")
        path = Path(value).expanduser()
        result[key] = (path if path.is_absolute() else config.parent / path).resolve()
    return result
