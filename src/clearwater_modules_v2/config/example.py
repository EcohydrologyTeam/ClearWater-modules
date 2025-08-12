"""Write out an example configuration file"""

import yaml
from pathlib import Path

config = {
    "model": {
        "start_time": "2022-05-14 00:00:00",
        "end_time": "2022-05-15 00:00:00",
        "time_step": "30s",
    },
    "processes": [
        {"riverine": {"configuration_path": "riverine_config.yml"}},
        {"temperature": {"wind_a": 1.3, "wind_b": 1.5, "wind_c": 3.0}},
    ],
    "sources": {
        "source_name_1": {
            "provider": "csv",
            "data": {"file_path": "path_to_csv"},
        }
    },
}


def save_example_config(path: Path):
    with open(path, "w") as f:
        yaml.dump(config, f)


if __name__ == "__main__":
    save_example_config(Path("example_config.yml"))
