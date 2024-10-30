import sys
import os
import random
import numpy as np
import argparse
import logging
import yaml
from outputfunc import *
import ast
from test_module import *

sys.path.append(".")

def parse_arguments():
    """
        Parse command-line arguments.

        Returns:
            argparse.Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--config_file", required=True, help="The path to the config file in yaml format. ")
    parser.add_argument("--inpdir", type=str, help="input file path")
    # parser.add_argument("--library", type=str, help="path for library")
    parser.add_argument("--device", choices=["cuda", "cpu"], default="cpu")
    parser.add_argument("--logdir", required=True, type=str, help="info.log will be generated")
    # parser.add_argument("--seed", type=int, help="random seed for weight initialization")
    parser.add_argument("--output_file_name", required=True, type=str, help="output file name", default="pydoseia.txt")
    directory = parser.parse_args().logdir
    if not os.path.exists(directory):
        os.makedirs(directory)
    return parser.parse_args()


def parse_config(configfile: str):
    """
        Parse configuration file.

        Args:
            configfile (str): Path to the configuration file in YAML format.

        Returns:
            dict: Parsed configuration data.
    """
    with open(configfile, "r") as f:
        data = yaml.load(f, Loader=yaml.FullLoader)
        # yaml.load(f, Loader=yaml.FullLoader)
        return convert_none_to_str(data)


def convert_none_to_str(data):
    """
        Convert None values to string 'None' in a nested data structure.

        Args:
            data: Data structure to process.

        Returns:
            str or dict or list: Processed data structure.
    """
    if isinstance(data, list):
        data[:] = [convert_none_to_str(i) for i in data]
    elif isinstance(data, dict):
        for k, v in data.items():
            data[k] = convert_none_to_str(v)
    return 'None' if data is None else data


def main():
    """
        Main function to parse arguments, configure logging, parse configuration file, perform dose calculation, and generate output.
    """
    args = parse_arguments()
    # seed_torch(seed=int(args.seed))
    log_file_name = args.config_file.split('.yaml')[0]

    # Configure logging in the main script
    logging.basicConfig(
        filename=f'{log_file_name}_info.log',
        level=logging.INFO,
        filemode="w+",
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

    # Add a handler to log to stdout
    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setLevel(logging.INFO)
    logging.getLogger().addHandler(stdout_handler)

    # logging.basicConfig(filename='%s_info.log' % log_file_name, level=logging.INFO, filemode="w+")
    # logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))
    # logging.getLogger("main").info(f"cuda available: {torch.cuda.is_available()}")
    config = parse_config(args.config_file)
    # metfunc = MetFunc(config=config, device=args.device, logdir=args.logdir)
    output_func = OutputFunc(config=config, device=args.device, log_file_name=log_file_name, logdir=args.logdir)
    # all_results = output_func.dose_calculation_script()

    if not config['run_dose_computation']:
        dilution_factor_sectorwise_all_distances = output_func.dose_calculation_script()
        output_func.output_to_txt(dilution_factor_sectorwise_all_distances, filename=args.output_file_name,
                                  DCFs=None, DOSES=None, INGESTION_DOSES=None, PLUME_DOSES=None)

    if config['run_dose_computation']:
        if config['run_plume_shine_dose']:
            DCFs, dilution_factor_sectorwise_all_distances, DOSES, INGESTION_DOSES, PLUME_DOSES = output_func.dose_calculation_script()
            output_func.output_to_txt(dilution_factor_sectorwise_all_distances, filename=args.output_file_name, DCFs=DCFs,
                                  DOSES=DOSES, INGESTION_DOSES=INGESTION_DOSES, PLUME_DOSES=PLUME_DOSES)
        else:
            DCFs, dilution_factor_sectorwise_all_distances, DOSES, INGESTION_DOSES = output_func.dose_calculation_script()
            output_func.output_to_txt(dilution_factor_sectorwise_all_distances, filename=args.output_file_name,
                                  DCFs=DCFs, DOSES=DOSES, INGESTION_DOSES=INGESTION_DOSES, PLUME_DOSES=None)

    logging.getLogger("main").info("path of input file: {config_desc}".format(config_desc=args.config_file))
    logging.getLogger("main").info(
        "output file name: {output_file_name}".format(output_file_name=args.output_file_name))


if __name__ == "__main__":
    main()
