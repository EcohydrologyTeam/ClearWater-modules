"""Reads a model configuration file"""

import yaml


def read_config(config_path: str) -> dict:
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)
    return config
