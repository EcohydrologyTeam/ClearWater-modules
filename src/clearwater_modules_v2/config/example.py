"""Write out an example configuration file"""

from pathlib import Path
import sys
import shutil


def save_example_config(target_path: Path) -> None:
    source_path = Path(__file__).parent / "example_config.yml"
    shutil.copy(source_path, target_path)


def main(*args) -> None:
    # if no directory provided, use current directory
    if len(args) < 2:
        return save_example_config(Path.cwd() / "example_config.yml")
    path = Path(args[1])
    # if user provided a file, confirm that it is a yaml file
    if path.suffix == ".yml" or path.suffix == ".yaml":
        return save_example_config(path)
    raise ValueError(
        "Provided pathway must be a valid `.yaml` or `yml` file or a directory location."
    )


if __name__ == "__main__":
    main(*sys.argv[1:])
