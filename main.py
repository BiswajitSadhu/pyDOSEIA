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
import _pickle as cPickle
import pandas as pd

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


def process_plume_doses(loaded_data, radionuclides, distances=[100, 200]):

    """
    Processes plume dose data from a pickle file and returns a structured DataFrame.

    Parameters:
    - pickle_file (str): Path to the pickle file containing plume dose data.
    - radionuclides (list): List of radionuclide names.
    - distances (list, optional): List of distance labels. Default is [100, 200].

    Returns:
    - pd.DataFrame: Processed plume dose DataFrame.
    """
    # Load the data
    # with open(pickle_file, "rb") as f:
    #    loaded_data = cPickle.load(f)

    PDS = loaded_data['PLUME_DOSES']
    sectors = ['N', 'NNE', 'NE', 'ENE', 'E', 'ESE', 'SE',
               'SSE', 'S', 'SSW', 'SW', 'WSW', 'W', 'WNW', 'NW', 'NNW']

    num_radionuclides = len(radionuclides)

    # Initialize lists to store data
    psd = []
    distance_labels = []
    radionuclide_labels = []

    for dist_idx, each in enumerate(PDS):
        each = np.array(each, dtype=object).sum(axis=0)  # Sum along axis 0
        each_flat = [x.ravel() for x in each]  # Flatten each sub-array
        each_concat = np.concatenate(each_flat)  # Concatenate all flattened values
        each = each_concat.reshape(-1, 16)  # Reshape to have 16 sectors per row

        # Append reshaped data
        psd.append(each)

        # Assign distance labels
        distance_labels.extend([distances[dist_idx % len(distances)]] * num_radionuclides)

        # Assign Radionuclide labels
        radionuclide_labels.extend(radionuclides)

    # Stack all `each` arrays row-wise
    psd_stacked = np.vstack(psd)

    # Convert to DataFrame
    df_plume_dose = pd.DataFrame(psd_stacked, dtype=float, columns=sectors)

    # Add Distance and Radionuclide columns
    df_plume_dose.insert(0, "Distance (m)", distance_labels)
    df_plume_dose.insert(1, "Radionuclide", radionuclide_labels)

    # Compute max values and sectors
    df_plume_dose['Maximum'] = df_plume_dose.iloc[:, 2:].max(axis=1)
    df_plume_dose['Maximum Sector'] = df_plume_dose.iloc[:, 2:].idxmax(axis=1)

    # save to .csv
    df_plume_dose.to_csv('summary_plume_dose.csv')

    return df_plume_dose


def reshape_ingestion_dose_data(ingestion_doses, ages, distances, radionuclides, plant_boundary_dist):
    """
    Reshapes a 4D ingestion dose array (with three ingestion routes) into a structured DataFrame.

    Parameters:
    - ingestion_doses: np.array of shape (num_ages, num_distances, num_radionuclides, 3).
    - ages: List of age values.
    - distances: List of distance values.
    - radionuclides: List of radionuclide names.

    Returns:
    - A pandas DataFrame with Age, Distance, Radionuclide, Veg Dose, Milk Dose, Meat Dose, and Total Dose.
      Includes summed rows per distance.
    """

    if plant_boundary_dist not in distances:
        distances = distances + [plant_boundary_dist]

    # data_rows = []

    # Extract dimensions
    num_distances, num_ages, num_radionuclides, num_routes = ingestion_doses.shape
    # print("SHAPE_INGESTION:::", ingestion_doses.shape)
    # Validate input dimensions
    assert len(ages) == num_ages, "Mismatch in age list length and array dimensions"
    assert len(distances) == num_distances, "Mismatch in distance list length and array dimensions"
    assert len(radionuclides) == num_radionuclides, "Mismatch in radionuclide list length and array dimensions"
    assert num_routes == 3, "Ingestion dose array should have 3 ingestion routes (Veg, Milk, Meat)"

    df_list = []  # Store dataframes before concatenation

    # Iterate over all dimensions
    for i, dist in enumerate(distances):  # Iterate over ages first
        for j, ag in enumerate(ages):  # Iterate over distances
            dose_entries = []  # Store individual rows for each distance

            for k, radionuclide in enumerate(radionuclides):  # Iterate over radionuclides
                veg_dose = ingestion_doses[i, j, k, 0]  # Dose from vegetables
                milk_dose = ingestion_doses[i, j, k, 1]  # Dose from milk
                meat_dose = ingestion_doses[i, j, k, 2]  # Dose from meat
                total_dose = veg_dose + milk_dose + meat_dose  # Total ingestion dose

                dose_entries.append([dist, ag, radionuclide, veg_dose, milk_dose, meat_dose, total_dose])

            # Convert to DataFrame
            df_temp = pd.DataFrame(dose_entries, columns=['Distance (m)', 'Age (y)', 'Radionuclide',
                                                          'Veg Dose', 'Milk Dose', 'Meat Dose', 'Total'])

            # Compute summed row
            sum_row = df_temp[['Veg Dose', 'Milk Dose', 'Meat Dose', 'Total']].sum()
            sum_row = pd.DataFrame([{
                'Distance (m)': dist,
                'Age (y)': ag,
                'Radionuclide': 'SUM',  # Marker for summed row
                'Veg Dose': sum_row['Veg Dose'],
                'Milk Dose': sum_row['Milk Dose'],
                'Meat Dose': sum_row['Meat Dose'],
                'Total': sum_row['Total']
            }])

            # Append summed row
            df_temp = pd.concat([df_temp, sum_row], ignore_index=True)

            # Store the results
            df_list.append(df_temp)

    # Combine all DataFrames
    df = pd.concat(df_list, ignore_index=True)

    # Sort the DataFrame by Age to ensure Age 1 comes before Age 18
    df = df.sort_values(by=['Age (y)'], ascending=True).reset_index(drop=True)

    # Convert scientific notation for readability
    pd.options.display.float_format = '{:.3e}'.format

    return df

# Example Usage
# ages = [1, 18]  # Age groups
# distances = [100, 200]  # Distance values
# radionuclides = ['Sr-90', 'Cs-137', 'H-3']  # Radionuclide types


def reshape_dose_data(doses_array, df_ing, ages, distances, radionuclides, plant_boundary_dist):
    """
    Reshapes a 4D dose array into a structured DataFrame and adds summed rows per distance.
    Also adds ingestion dose from df_ing.

    Parameters:
    - doses_array: np.array of shape (num_ages, num_distances, dose_types, num_radionuclides).
    - df_ing: DataFrame containing ingestion doses.
    - ages: List of age values.
    - distances: List of distance values.
    - radionuclides: List of radionuclide names.
    - plant_boundary_dist: distance of plant boundary (metre)

    Returns:
    - df: DataFrame with individual radionuclide doses.
    - df_sum_only: DataFrame with summed doses per distance, including ingestion dose.
    """

    if plant_boundary_dist not in distances:
        distances = distances + [plant_boundary_dist]

    df_list = []  # Store dataframes before concatenation
    sum_list = []  # Store summed rows

    # Extract dimensions
    num_distances, num_ages, dose_types, num_radionuclides = doses_array.shape

    # Validate input dimensions
    assert len(ages) == num_ages, "Mismatch in age list length and array dimensions"
    assert len(distances) == num_distances, "Mismatch in distance list length and array dimensions"
    assert len(radionuclides) == num_radionuclides, "Mismatch in radionuclide list length and array dimensions"
    assert dose_types == 3, "Dose array should have 3 dose types (Inhalation, Ground-shine, Submersion)"
    
    # Iterate over all dimensions
    for i, dist in enumerate(distances):
        for j, ag in enumerate(ages):
            dose_entries = []  # Store individual rows for each distance

            for k, radionuclide in enumerate(radionuclides):  # Iterate over radionuclides
                inhalation = doses_array[i, j, 0, k]  
                ground_shine = doses_array[i, j, 1, k]  
                submersion = doses_array[i, j, 2, k]  
                total_dose = inhalation + ground_shine + submersion

                dose_entries.append([dist, ag, radionuclide, inhalation, ground_shine, submersion, total_dose])

            # Convert to DataFrame
            df_temp = pd.DataFrame(dose_entries, columns=['Age (y)', 'Distance (m)', 'Radionuclide', 
                                                          'Inhalation dose (mSv)', 'Ground-shine (mSv)', 
                                                          'Submersion (mSv)', 'Total (mSv)'])
            
            # Compute summed row (without Radionuclide column, and removing Total)
            sum_row = df_temp[['Inhalation dose (mSv)', 'Ground-shine (mSv)', 'Submersion (mSv)']].sum()
            sum_row = {
                'Age (y)': ag,
                'Distance (m)': dist,
                'Inhalation dose (mSv)': sum_row['Inhalation dose (mSv)'], 
                'Ground-shine (mSv)': sum_row['Ground-shine (mSv)'], 
                'Submersion (mSv)': sum_row['Submersion (mSv)']
            }

            # Append summed row separately
            sum_list.append(sum_row)

            # Store the results
            df_list.append(df_temp)

    # Combine all DataFrames
    df = pd.concat(df_list, ignore_index=True)
    df_sum_only = pd.DataFrame(sum_list)  # DataFrame with only SUM rows (without Radionuclide and Total)

    # Add ingestion dose from df_ing
    df_ing_total = df_ing.groupby(['Age (y)', 'Distance (m)'])['Total'].sum().reset_index()
    df_ing_total.rename(columns={'Total': 'Ingestion dose (mSv)'}, inplace=True)

    # Merge ingestion dose into df_sum_only
    df_sum_only = df_sum_only.merge(df_ing_total, on=['Age (y)', 'Distance (m)'], how='left')

    # Compute final total including ingestion
    df_sum_only['Total Dose (mSv)'] = df_sum_only[['Inhalation dose (mSv)', 'Ground-shine (mSv)', 'Submersion (mSv)', 'Ingestion dose (mSv)']].sum(axis=1)

    # Sort the DataFrames
    df = df.sort_values(by=['Age (y)']).reset_index(drop=True)
    df_sum_only = df_sum_only.sort_values(by=['Age (y)']).reset_index(drop=True)

    # Convert scientific notation for readability
    pd.options.display.float_format = '{:.3e}'.format
    
    return df, df_sum_only

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
        ages = config['age_group']
        plant_boundary_dist = config['plant_boundary']
        distances = config['downwind_distances']
        radionuclides = config['rads_list']

        if config['run_plume_shine_dose']:
            DCFs, dilution_factor_sectorwise_all_distances, DOSES, INGESTION_DOSES, PLUME_DOSES = output_func.dose_calculation_script()
            output_func.output_to_txt(dilution_factor_sectorwise_all_distances, filename=args.output_file_name, DCFs=DCFs,
                                  DOSES=DOSES, INGESTION_DOSES=INGESTION_DOSES, PLUME_DOSES=PLUME_DOSES)
        
            # Save multiple arrays
            data = {"DOSES": DOSES, "ING_DOSES": INGESTION_DOSES, "PLUME_DOSES": PLUME_DOSES}
            
            with open("doses_ing_doses_ps_doses.pkl", "wb") as f:
                cPickle.dump(data, f)



            # INGESTION DOSE (three routes)
            df_ing = reshape_ingestion_dose_data(data['ING_DOSES'], ages, distances, radionuclides, plant_boundary_dist)
            df_ing.to_csv('summary_detailed_ingestion_dose.csv')

            # INH_GS_SUBS DOSES
            df, df_sum_only = reshape_dose_data(data['DOSES'], df_ing, ages, distances, radionuclides, plant_boundary_dist)
            
            df.to_csv('summary_detailed_inh_gs_sub_dose.csv')
            df_sum_only.to_csv('summary_summed_inh_gs_sub_dose.csv')

            print("Full DataFrame:\n", df)
            print("\nSummed DataFrame (with Ingestion dose and Total Dose):\n", df_sum_only)
            
            # PLUME DOSE
            if config['long_term_release']:
                df_plume = process_plume_doses(data, radionuclides, distances)
                df_plume.to_csv('summary_detailed_plume_dose.csv')    

        else:
            DCFs, dilution_factor_sectorwise_all_distances, DOSES, INGESTION_DOSES = output_func.dose_calculation_script()
            output_func.output_to_txt(dilution_factor_sectorwise_all_distances, filename=args.output_file_name,
                                  DCFs=DCFs, DOSES=DOSES, INGESTION_DOSES=INGESTION_DOSES, PLUME_DOSES=None)

            # Save multiple arrays
            data = {"DOSES": DOSES, "ING_DOSES": INGESTION_DOSES}

            with open("doses_ing_doses.pkl", "wb") as f:
                cPickle.dump(data, f)

            # INGESTION DOSE (three routes)
            df_ing = reshape_ingestion_dose_data(data['ING_DOSES'], ages, distances, radionuclides, plant_boundary_dist)
            df_ing.to_csv('summary_detailed_ingestion_dose.csv')

            # INH_GS_SUBS DOSES
            df, df_sum_only = reshape_dose_data(data['DOSES'], df_ing, ages, distances, radionuclides, plant_boundary_dist)

            df.to_csv('summary_detailed_inh_gs_sub_dose.csv')
            df_sum_only.to_csv('summary_summed_inh_gs_sub_dose.csv')

            print("Full DataFrame:\n", df)
            print("\nSummed DataFrame (with Ingestion dose and Total Dose):\n", df_sum_only)

            with open("doses_ing_doses.pkl", "wb") as f:
                cPickle.dump(data, f)

    logging.getLogger("main").info("path of input file: {config_desc}".format(config_desc=args.config_file))
    logging.getLogger("main").info(
        "output file name: {output_file_name}".format(output_file_name=args.output_file_name))


if __name__ == "__main__":
    main()
