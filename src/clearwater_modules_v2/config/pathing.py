from pathlib import Path, PureWindowsPath
from typing import Optional

# TODO: move to ClearWater-data, because this is shared with riverine

def resolve_path(path: str|Path, repo_path: Optional[str|Path] = None):
    """Resolves filepath from configuration file"""
    # convert windows path string, if supplied
    path = Path(PureWindowsPath(path))

    if path.is_absolute():
        absolute_path = path
    else:
        if repo_path is None:
            # Repo path, relative to this module
            repo_path = Path(__file__).parent.parent.parent.parent
        absolute_path = repo_path / path

    validate_path(absolute_path)
    return absolute_path
    

def validate_path(path: Path):
    """Validate if path exists"""
    if not path.exists():
        raise FileNotFoundError(f"Path does not exist: {path}")
    