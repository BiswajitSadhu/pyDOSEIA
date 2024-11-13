import ast
import yaml
import os.path
from colorama import init, Fore, Back, Style
import pandas as pd
from auto_input_generator_funcs_class import *
import subprocess
import random

# must have xlrd>=2.0.1 to read excel file
init()
color = Fore.GREEN
color_error = Fore.RED

InpGenFunc = InpGenFunc()


def auto_input_generator():
    run_dose, input_dict = InpGenFunc.get_run_dose_computation()
    # DOSE + DILUTION FACTOR
    if run_dose == 'Y':
        # DOSE RELATED
        rads_list, input_dict = InpGenFunc.get_rads_list()
        print('\n')
        type_rad, input_dict = InpGenFunc.get_type_rad(rads_list)
        print('\n')
        get_consider_progeny = InpGenFunc.get_consider_progeny()
        print('\n')
        # only for ground shine dose computation
        input_dict = InpGenFunc.get_weathering_corr_gs()
        # print('\n')
        # only for ground shine dose computation
        # input_dict = InpGenFunc.get_weathering_corr_ingestion() ; CHECK later
        print('\n')
        exposure_period, input_dict = InpGenFunc.get_exposure_period_for_ground_shine()
        print('\n')
        input_dict = InpGenFunc.get_run_plume_shine_dose()
        print('\n')
        # age list for calculation; must be a list e.g. [1,18]; 1 for infant and 18 for adult.
        age_group, input_dict = InpGenFunc.get_age_group()
        print('\n')
        ad = len([_ for _ in age_group if _ >= 18])
        if 1 in age_group and ad > 0:
            inges_param_dict_adult, input_dict = InpGenFunc.get_inges_param_dict_adult()
            print('\n')
            inges_param_dict_infant, input_dict = InpGenFunc.get_inges_param_dict_infant()
        if 1 in age_group and ad == 0:
            inges_param_dict_infant, input_dict = InpGenFunc.get_inges_param_dict_infant()
        if 1 not in age_group and ad > 0:
            inges_param_dict_adult, input_dict = InpGenFunc.get_inges_param_dict_adult()

        # DILUTION FACTOR COMPUTATION
        print('\n')
        release_height, input_dict = InpGenFunc.get_release_height()
        print('\n')
        downwind_distances, input_dict = InpGenFunc.get_downwind_distances()
        print('\n')
        plant_boundary, input_dict = InpGenFunc.get_plant_boundary(downwind_distances)
        print('\n')
        # Type of release scenario
        release_scenario, input_dict = InpGenFunc.get_release_scenario()
        print('\n')
        if release_scenario == 1:
            print("you opted for Long Term Continuous Release")
            # ask user about Y and Z values for using master equation as per need.
            print('\n')
            input_dict = InpGenFunc.get_master_eq_continuous_plume()
            print('\n')
            annual_discharge_bq_rad_list, input_dict = InpGenFunc.get_annual_discharge_bq_rad_list(rads_list)
            print('\n')
            have_dilution_factor, input_dict = InpGenFunc.get_details_have_dilution_factor()
            print('\n')
            # measurement height of meteorological data
            measurement_height, input_dict = InpGenFunc.get_measurement_height()
            if have_dilution_factor == 'N':
                have_met_data_bool, input_dict = InpGenFunc.get_have_met_data_release_scenario_1(release_scenario)
                print('have_met_data_bool:', have_met_data_bool)
                if have_met_data_bool:
                    print('\n')
                    path_met_file, input_dict = InpGenFunc.get_path_met_file()
                    print('\n')
                    operation_time, input_dict = InpGenFunc.get_operation_time()
                    print('\n')
                    excel_sheet_name, input_dict = InpGenFunc.get_excel_sheet_name()
                    print('\n')
                    column_names, input_dict = InpGenFunc.get_column_names()
                    print('\n')
                    num_days, input_dict = InpGenFunc.get_num_days(excel_sheet_name)
                    print('\n')
                    # calm correction (Yes/No)
                    calm_correction, input_dict = InpGenFunc.get_calm_correction()
                    print('\n')
                # no met data but long-term release; use conservative approach
                else:
                    print('\n')
                    print("# no met data available but long-term release; using conservative approach")
                    like_to_scale_with_mean_speed, input_dict = InpGenFunc.scale_dilution_factor_with_user_defined_mean_speed()
                    if like_to_scale_with_mean_speed == 'Y':
                        ask_mean_speed_data, input_dict = InpGenFunc.get_mean_speed_stab_cat_wise()
                    elif like_to_scale_with_mean_speed == 'N':
                        print(
                            "dilution factor won't be scaled with meteorological data as the data is not provided by the "
                            "user. Mean wind speed = 1 will be used.")
                    print('\n')


            elif have_dilution_factor == 'Y':
                input_dict = InpGenFunc.get_have_met_data(have_dilution_factor, release_scenario)
                print('\n')
                max_dilution_factor, input_dict = InpGenFunc.get_max_dilution_factor(downwind_distances)
            else:
                print("input must be Y or N.")

        # RELEASE SCENARIO 2; INSTANTENOUS RELEASE
        elif release_scenario == 2:
            # print("you opted for Instantenous Release ")
            # consumption_time_food, input_dict = InpGenFunc.get_consumption_time_food()
            print('\n')
            instantaneous_release_bq_list, input_dict = InpGenFunc.get_instantaneous_release_bq_list(rads_list)
            print('\n')
            have_dilution_factor, input_dict = InpGenFunc.get_details_have_dilution_factor()
            print('\n')
            if have_dilution_factor == 'N':
                print('\n')
                # ask user about Y and Z values for using master equation as per need.
                input_dict = InpGenFunc.get_master_eq_single_plume()
                have_met_data, input_dict = InpGenFunc.get_have_met_data_release_scenario_2(release_scenario)
                if have_met_data == 'Y':
                    print('\n')
                    operation_time, input_dict = InpGenFunc.get_operation_time()
                    print('\n')
                    path_met_file, input_dict = InpGenFunc.get_path_met_file()
                    print('\n')
                    operation_time, input_dict = InpGenFunc.get_operation_time()
                    print('\n')
                    excel_sheet_name, input_dict = InpGenFunc.get_excel_sheet_name()
                    print('\n')
                    column_names, input_dict = InpGenFunc.get_column_names()
                    print('\n')
                    num_days, input_dict = InpGenFunc.get_num_days(excel_sheet_name)
                    print('\n')
                    # measurement height of meteorological data
                    measurement_height, input_dict = InpGenFunc.get_measurement_height()
                    print('\n')
                    # calm correction (Yes/No)
                    print('\n')
                    calm_correction, input_dict = InpGenFunc.get_calm_correction()

                # no precomputed dilution factor and no metereological data: dilution factor
                # will be computed but won't be scaled
                elif have_met_data == 'N':
                    print('\n')
                    # measurement height of meteorological data; N:B Wind speeds
                    # measured at the 10 m level should be used for those times when the
                    # effluent plume is considered to be a ground level release. (ref: IAEA TECDOC 379, page 60)
                    measurement_height, input_dict = InpGenFunc.get_measurement_height()
                    print('\n')
                    like_to_scale_with_mean_speed, input_dict = InpGenFunc.scale_dilution_factor_with_user_defined_mean_speed()
                    if like_to_scale_with_mean_speed == 'Y':
                        ask_mean_speed_data, input_dict = InpGenFunc.get_mean_speed_stab_cat_wise()
                    elif like_to_scale_with_mean_speed == 'N':
                        print(
                            "dilution factor won't be scaled with meteorological data as the data is not provided by the "
                            "user. Mean wind speed = 1 will be used.")

                    # awesome, the precomputed dilution factor will be used
            elif have_dilution_factor == 'Y':
                print('\n')
                input_dict = InpGenFunc.get_have_met_data(have_dilution_factor, release_scenario)
                print('\n')
                max_dilution_factor, input_dict = InpGenFunc.get_max_dilution_factor(downwind_distances)
        else:
            # raise Exception('Please press either 1 or 2.')
            print("Please press either 1 or 2.")

        if 'C-14' in rads_list or 'H-3' in rads_list:
            print('\n')
            veg_type_list, input_dict = InpGenFunc.get_veg_type_list_h3c14()
            # options: arctic, Mediterranean, meritime, Continental
            print('\n')
            climate, input_dict = InpGenFunc.get_climate_h3c14()
            print('\n')
            animal_feed_type = InpGenFunc.get_animal_feed_type_h3c14()
            print('\n')
            if 'H-3' in rads_list:
                print('\n')
                animal_product_list_for_tritium, input_dict = InpGenFunc.get_animal_product_list_for_tritium()
            if 'C-14' in rads_list:
                print('\n')
                animal_product_list_for_C14, input_dict = InpGenFunc.get_animal_product_list_for_C14()
        print('\n')
        soiltype, input_dict = InpGenFunc.get_soiltype()
        print('\n')
        inges_param_dict, input_dict = InpGenFunc.get_inges_param_dict()

    # ONLY DILUTION FACTOR
    if run_dose == 'N':
        # print('\n')
        # DILUTION FACTOR COMPUTATION
        # rads_list, input_dict = InpGenFunc.get_rads_list()
        print('\n')
        release_height, input_dict = InpGenFunc.get_release_height()
        print('\n')
        downwind_distances, input_dict = InpGenFunc.get_downwind_distances()
        print('\n')
        plant_boundary, input_dict = InpGenFunc.get_plant_boundary(downwind_distances)
        print('\n')
        # Type of release scenario
        release_scenario, input_dict = InpGenFunc.get_release_scenario()
        print('\n')
        # measurement height of meteorological data
        measurement_height, input_dict = InpGenFunc.get_measurement_height()
        print('\n')
        if release_scenario == 1:
            # print("you opted for Long Term Continuous Release")
            # ask user about Y and Z values for using master equation as per need.
            print('\n')
            input_dict = InpGenFunc.get_master_eq_continuous_plume()
            # print('\n')
            # annual_discharge_bq_rad_list, input_dict = InpGenFunc.get_annual_discharge_bq_rad_list(rads_list)
            print('\n')
            have_dilution_factor, input_dict = InpGenFunc.get_details_have_dilution_factor()
            print('\n')
            if have_dilution_factor == 'N':
                have_met_data_bool, input_dict = InpGenFunc.get_have_met_data_release_scenario_1(release_scenario)
                print('have_met_data_bool:', have_met_data_bool)
                if have_met_data_bool:
                    print('\n')
                    path_met_file, input_dict = InpGenFunc.get_path_met_file()
                    print('\n')
                    operation_time, input_dict = InpGenFunc.get_operation_time()
                    print('\n')
                    excel_sheet_name, input_dict = InpGenFunc.get_excel_sheet_name()
                    print('\n')
                    column_names, input_dict = InpGenFunc.get_column_names()
                    print('\n')
                    num_days, input_dict = InpGenFunc.get_num_days(excel_sheet_name)
                    print('\n')
                    # calm correction (Yes/No)
                    calm_correction, input_dict = InpGenFunc.get_calm_correction()
                    print('\n')
                # no met data but long-term release; use conservative approach
                else:
                    print('\n')
                    print("# no met data available but long-term release; using conservative approach")
                    like_to_scale_with_mean_speed, input_dict = InpGenFunc.scale_dilution_factor_with_user_defined_mean_speed()
                    if like_to_scale_with_mean_speed == 'Y':
                        ask_mean_speed_data, input_dict = InpGenFunc.get_mean_speed_stab_cat_wise()
                    elif like_to_scale_with_mean_speed == 'N':
                        print(
                            "dilution factor won't be scaled with meteorological data as the data is not provided by the "
                            "user. Mean wind speed = 1 will be used.")
                    print('\n')


            elif have_dilution_factor == 'Y':
                input_dict = InpGenFunc.get_have_met_data(have_dilution_factor, release_scenario)
                print('\n')
                max_dilution_factor, input_dict = InpGenFunc.get_max_dilution_factor(downwind_distances)
            else:
                print("input must be Y or N.")

        # RELEASE SCENARIO 2; INSTANTENOUS RELEASE
        elif release_scenario == 2:
            # print("you opted for Instantenous Release ")
            # consumption_time_food, input_dict = InpGenFunc.get_consumption_time_food()
            # print('\n')
            # instantaneous_release_bq_list, input_dict = InpGenFunc.get_instantaneous_release_bq_list(rads_list)
            print('\n')
            have_dilution_factor, input_dict = InpGenFunc.get_details_have_dilution_factor()
            print('\n')
            if have_dilution_factor == 'N':
                print('\n')
                # ask user about Y and Z values for using master equation as per need.
                input_dict = InpGenFunc.get_master_eq_single_plume()
                have_met_data, input_dict = InpGenFunc.get_have_met_data_release_scenario_2(release_scenario)
                if have_met_data == 'Y':
                    print('\n')
                    operation_time, input_dict = InpGenFunc.get_operation_time()
                    print('\n')
                    path_met_file, input_dict = InpGenFunc.get_path_met_file()
                    print('\n')
                    operation_time, input_dict = InpGenFunc.get_operation_time()
                    print('\n')
                    excel_sheet_name, input_dict = InpGenFunc.get_excel_sheet_name()
                    print('\n')
                    column_names, input_dict = InpGenFunc.get_column_names()
                    print('\n')
                    num_days, input_dict = InpGenFunc.get_num_days(excel_sheet_name)
                    print('\n')
                    # measurement height of meteorological data
                    # measurement_height, input_dict = InpGenFunc.get_measurement_height()
                    # print('\n')
                    # calm correction (Yes/No)
                    print('\n')
                    calm_correction, input_dict = InpGenFunc.get_calm_correction()

                # no precomputed dilution factor and no metereological data: dilution factor
                # will be computed but won't be scaled
                elif have_met_data == 'N':
                    print('\n')
                    # measurement height of meteorological data; N:B Wind speeds
                    # measured at the 10 m level should be used for those times when the
                    # effluent plume is considered to be a ground level release. (ref: IAEA TECDOC 379, page 60)
                    # measurement_height, input_dict = InpGenFunc.get_measurement_height()
                    # print('\n')
                    like_to_scale_with_mean_speed, input_dict = InpGenFunc.scale_dilution_factor_with_user_defined_mean_speed()
                    if like_to_scale_with_mean_speed == 'Y':
                        ask_mean_speed_data, input_dict = InpGenFunc.get_mean_speed_stab_cat_wise()
                    elif like_to_scale_with_mean_speed == 'N':
                        print(
                            "dilution factor won't be scaled with meteorological data as the data is not provided by the "
                            "user. Mean wind speed = 1 will be used.")

                    # awesome, the precomputed dilution factor will be used
            elif have_dilution_factor == 'Y':
                print('\n')
                input_dict = InpGenFunc.get_have_met_data(have_dilution_factor, release_scenario)
                print('\n')
                max_dilution_factor, input_dict = InpGenFunc.get_max_dilution_factor(downwind_distances)
        else:
            # raise Exception('Please press either 1 or 2.')
            print("Please press either 1 or 2.")

    # generating INPUT FILE, default = False
    pickle_it = InpGenFunc.get_pickle_it()
    # set run_ps_dose_parallel; default = True
    run_ps_dose_parallel = InpGenFunc.get_run_ps_dose_parallel()
    if run_ps_dose_parallel:
        print("PLUME SHINE DOSE WILL BE CALCULATED THROUGH PARALLEL PROCESSORS..")
    else:
        print("PLUME SHINE DOSE WILL BE CALCULATED THROUGH SINGLE PROCESSOR..")

    with open('%s/%s.yaml' % (logdir_name, input_file_name), 'w') as infile:
        yaml.safe_dump(InpGenFunc.input_dict, infile, default_flow_style=True)
    print("generated content of input file:", InpGenFunc.input_dict)



if __name__ == "__main__":
    logdir_name, input_dict = InpGenFunc.get_logdir_name()

    if not os.path.exists('%s' % logdir_name):
        os.makedirs('%s' % logdir_name)
    input_file_name, file_exists, input_dict = InpGenFunc.get_input_file_name(logdir_name)

    output_file_name, input_dict = InpGenFunc.get_output_file_name()

    # subprocess.call("mv %s.yaml %s/" % (input_file_name, input_dict['logdir_name']), shell=True)
    if not file_exists:
        auto_input_generator()

    # check input file for consistency in entries
    InpGenFunc.check_existing_input_entries()

    print("_________________________________________________________________________________")
    print("That's all we needed. Now sit back and relax while we perform the computation :).")

    # Specify the log file path
    log_file_path = f"{input_dict['logdir_name']}/{input_dict['output_file_name']}.log"
    '''
    # Open the file in write mode and pass it to subprocess
    with open(log_file_path, "w") as log_file:
        subprocess.call(
            "python main.py --config_file %s/%s.yaml --output_file_name %s/%s --logdir %s" % (
                input_dict['logdir_name'],
                input_file_name,
                input_dict['logdir_name'],
                input_dict['output_file_name'],
                input_dict['logdir_name']
            ),
            shell=True,
            stdout=log_file,  # Redirect standard output to the file
            stderr=log_file  # Redirect standard error to the file
        )

    #subprocess.call("python main.py --config_file %s/%s.yaml --output_file_name %s/%s --logdir %s" % (
    #    input_dict['logdir_name'], input_file_name, input_dict['logdir_name'], input_dict['output_file_name'],
    #    input_dict['logdir_name']), shell=True)
    # Print completion message and quote on the screen
    print("\nCalculation finished successfully!")

    '''
    # Run the subprocess and capture output and errors
    command = (
        f"python main.py --config_file {input_dict['logdir_name']}/{input_file_name}.yaml "
        f"--output_file_name {input_dict['logdir_name']}/{input_dict['output_file_name']} "
        f"--logdir {input_dict['logdir_name']}"
    )

    # Execute the command and capture output
    result = subprocess.run(command, shell=True, capture_output=True, text=True)

    # Write output and error to the log file
    with open(log_file_path, "w") as log_file:
        log_file.write(result.stdout)
        log_file.write(result.stderr)

    # Check if the subprocess finished successfully
    if result.returncode == 0:
        print("\nCalculation finished successfully!")
    else:
        print("\nAn error occurred during the calculation. Please check the log file for more details.")
        print("Error details:\n", result.stderr)
        

    # Expanded list of inspiring quotes from Vedic literature
    vedic_quotes = [
        "Tat Tvam Asi (तत्त्वमसि) - Thou art that. - Chandogya Upanishad",
        "Aham Brahmasmi (अहं ब्रह्मास्मि) - I am Brahman. - Brihadaranyaka Upanishad",
        "Atman is Brahman. The self and the cosmic are one. - Advaita Vedanta",
        "Satyam eva jayate - Truth alone triumphs. - Mundaka Upanishad",
        "All that we are is the result of what we have thought. - Buddha",
        "He who sees all beings in his own self, and his own self in all beings, loses all fear. - Isa Upanishad",
        "There is no joy in the finite; there is joy only in the Infinite. - Chandogya Upanishad",
        "From the unreal, lead me to the real. From darkness, lead me to light. From death, lead me to immortality. - Brihadaranyaka Upanishad",
        "The self is not the body; the self is eternal, indestructible, and infinite. - Bhagavad Gita",
        "One who sees all beings as the Self and the Self in all beings has no hatred for any. - Isa Upanishad",
        "As the rivers flowing east and west merge in the sea and become one, so do all beings merge into the One. - Chandogya Upanishad",
        "He who knows himself to be the atman of all beings and all beings in the atman, hates none. - Isa Upanishad",
        "Where there is duality, there one sees another; but where everything has become one's own self, then who shall see whom? - Brihadaranyaka Upanishad",
        "Let your mind be clear and calm, like a serene lake. - Vedic Wisdom",
        "Be steadfast in performing your duty, O Arjuna, abandon attachment and remain balanced in success and failure. - Bhagavad Gita, Chapter 2, Verse 47",
        "The enlightened sees that Brahman is actionless and eternal, beyond good and evil. - Bhagavad Gita",
        "Knowing oneself to be Brahman, and seeing the self in all beings, one attains supreme peace. - Upanishadic Thought",
        "The wise see the same in all, be it a scholar, an outcast, an animal, or a sage. - Bhagavad Gita, Chapter 5, Verse 18",
        "Meditation brings wisdom; lack of meditation leaves ignorance. Know well what leads you forward and what holds you back. - Buddha",
        "Yoga is the journey of the self, through the self, to the self. - Bhagavad Gita",
        "The mind is restless, turbulent, powerful, and obstinate. It is like controlling the wind. - Bhagavad Gita, Chapter 6, Verse 34",
        "One who has control over the mind is tranquil in heat and cold, in pleasure and pain, and in honor and dishonor. - Bhagavad Gita, Chapter 6, Verse 7",
        "All creation is Brahman; and the universe is its play. - Upanishadic Thought",
        "The seeker of the Self should keep the mind calm, like the flame of a lamp in a windless place. - Bhagavad Gita, Chapter 6, Verse 19",
        "Look within. Be still. Free from fear and attachment, know the sweet joy of the way. - Dhammapada",
        "By pure thought, the wise man attains the supreme goal. - Vedic Wisdom"
    ]


    print("Vedas teach us:")
    print(random.choice(vedic_quotes))
