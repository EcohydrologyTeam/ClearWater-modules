from pathlib import Path, PureWindowsPath
from typing import Optional

# TODO: move to ClearWater-data, because this is shared with modules

def resolve_path(path: str|Path, project_path: Optional[str|Path] = None):
    """Resolves filepath from configuration file"""
    # convert windows path string, if supplied
    path = Path(PureWindowsPath(path))

    if path.is_absolute():
        absolute_path = path
    else:
        if project_path is None:
            # current working directory of the running process:
            # from wherever the user launched Python
            project_path = Path.cwd()
        absolute_path = project_path / path
    
    validate_path(absolute_path)
    return absolute_path
    

def validate_path(path: Path):
    """Validate if path exists"""
    if not path.exists():
        raise FileNotFoundError(f"Path does not exist: {path}")
    