import pandas as pd
from pathlib import Path
import yaml
from clearwater_modules_v2.config.init import init_from_file
from datetime import timedelta, datetime
import tempfile
import sys

# Coupled Riverine and Temperature for sumwere creek
MODULES_CONFIG_PATH = Path(r"C:\Users\sjordan\Github\ClearWater-modules\data_temp\sumwere_creek_fine_p49\modules_config.yml")
RIVERINE_CONFIG_PATH = Path(r"C:\Users\sjordan\Github\ClearWater-modules\data_temp\sumwere_creek_fine_p49\riverine_config.yml")
START_DATETIME = "2022-05-13 00:00:00"

def read_config(
    config_filepath: str | Path,
):
    with open(config_filepath, 'r') as file:
        model_config = yaml.safe_load(file)
    return model_config


def write_config(
    config: dict,
):
    tmp_file = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.yaml')
    yaml.safe_dump(config, tmp_file)
    tmp_file.close()
    return Path(tmp_file.name)


def model_updates(
        config: dict,
        chunk_size: timedelta,
        start_datetime: str | datetime,
        end_datetime: str | datetime,
    ):
    config['model']['start_datetime'] = start_datetime
    config['model']['end_datetime'] = end_datetime
    if chunk_size == "None":
        del config['model']['chunk_size']
    else:
        config['model']['chunk_size'] = chunk_size
    return config


def modify_config(
        chunk_size: timedelta,
        start_datetime: str | datetime,
        end_datetime: str | datetime,
    ):
    modules_config = read_config(MODULES_CONFIG_PATH)
    riverine_config = read_config(RIVERINE_CONFIG_PATH)

    riverine_config = model_updates(riverine_config, chunk_size, start_datetime, end_datetime)
    riverine_tmp_path = write_config(riverine_config)

    modules_config = model_updates(modules_config, chunk_size, start_datetime, end_datetime)
    modules_config['processes'][0]['riverine']['configuration_path'] = str(riverine_tmp_path)
    modules_tmp_path = write_config(modules_config)
    return modules_tmp_path


def run_model(modules_config: str | Path): 
    print("Initializing Model")
    model = init_from_file(modules_config)
    print("Running Model")
    model.run()


if __name__ == "__main__":

    import sys
    chunk_size, end_datetime = sys.argv[1], sys.argv[2]

    if chunk_size == "None":
        print("Running in non-chunked mode.")
    elif isinstance(chunk_size, str):
        try:
            pd.to_timedelta(chunk_size)
            print(f"Running in chunked mode with chunk size {chunk_size}")
        except(ValueError):
            raise ValueError(f"Invalid chunk size '{chunk_size}'. Must be a valid timedelta string.")
    else:
        raise TypeError("chunk_size must be None or a timedelta string")
    
    modules_tmp_path = modify_config(chunk_size, START_DATETIME, end_datetime)
    run_model(modules_tmp_path)

 