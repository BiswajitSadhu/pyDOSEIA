import ast
import yaml
import os.path
from colorama import init, Fore, Back, Style
import pandas as pd
import subprocess

init()
color = Fore.GREEN
color_error = Fore.RED


class InpGenFunc:
    def __init__(self):
        super().__init__()
        self.input_dict = {}

    def get_rads_list(self):
        rads_list = None
        while type(rads_list) != list:
            try:
                rads_list = ast.literal_eval(
                    input(color + "Provide the list of radionuclides (e.g. ['Co-60','H-3','C-14']): \n"))
                if type(rads_list) == list:
                    xls = pd.ExcelFile("library/Dose_ecerman_final.xlsx")
                    name = pd.read_excel(xls, "surface_dose")
                    name.dropna(axis=0, how='all', inplace=True)
                    for rad in rads_list:
                        search_string = '|'.join([rad])
                        df = name[name['Nuclide'].str.contains(search_string, na=False)]
                        df = df[df['Nuclide'].str.contains(r'(?:\s|^)%s(?:\s|$)' % search_string)]['Nuclide']
                        if df.empty:
                            print("could not idenfity the radionuclide {}".format(rad))
                            rads_list = None

                    print('Chosen set of radionuclide(s) is/are {}'.format(rads_list))
                else:
                    print(
                        color_error + "can't recognize input (not a list) (input formar should follow this: ['Co-60','H-3','C-14']).")
            # except Exception as ve:
            except Exception as ve:
                print(
                    color_error + "can't recognize input (not a list) (input formar should follow this: ['Co-60','H-3','C-14']).")

        self.input_dict['rads_list'] = rads_list
        return rads_list, self.input_dict

    def get_type_rad(self, rads_list):
        if type(rads_list) == list:
            # print("selected radioculides: {}".format(rads_list))
            element_list = [str(_).split('-')[0] for _ in rads_list]
            self.input_dict['element_list'] = element_list

            type_rad = None
            data_available = False
            while data_available == False or type(type_rad) != list or len(
                    [True for each in type_rad if each in ['Max', 'F', 'M', 'S']]) != len(rads_list):
                try:

                    type_rad = ast.literal_eval(input(
                        color + "Mention the type of radionculide to be used for fetching inhalation dose coeffecients. Format is: ['Max', 'F'] assuming two radionuclides are mentioned apriori. "
                                "\nFor the first radionuclide argument 'Max' ensures that the maximum value of inhalalation dose coefficients across all "
                                "the available types ('F', 'M', 'S') will be used for computation."
                                "While for the second type 'F' will be used. "
                                "\nPress 'D' for availing default option 'Max'. \n"))

                    if type_rad == 'D':
                        type_rad = ['Max' for each in rads_list]
                        data_available = True
                        self.input_dict['type_rad'] = type_rad
                        print("Your choice is: {}. It ensures that the maximum value of inhalalation dose coefficients"
                              " across all the available types ('F', 'M', 'S') will be used for computation.".format(
                            type_rad))


                    elif len(type_rad) != len(rads_list):
                        data_available = False
                        print(
                            "please provide input for all the radionuclides mentioned above. Format is: ['Max', 'F'] assuming two radionuclides are mentioned apriori. "
                            "\nFor the first radionuclide argument 'Max' ensures that the maximum value of inhalalation dose coefficients across all "
                            "the available types ('F', 'M', 'S') will be used for computation. While for the second type 'F' will be used ")

                    else:

                        print('Your choice is: {}'.format(type_rad))
                        master_file = 'library/RadioToxicityMaster.xls'
                        sheet_name = 'Inhalation CED Sv per Bq Public'
                        xls = pd.ExcelFile(master_file)
                        name = pd.read_excel(xls, sheet_name)
                        name.dropna(axis=0, how='all', inplace=True)

                        for ndx, rad in enumerate(rads_list):
                            search_string = '|'.join([rad])
                            df = name[name['Nuclide'] == search_string]
                            type_rad_elem = type_rad[ndx]
                            if type_rad_elem == 'Max':
                                # consider all types: F, M, S and take the max value as dcf
                                # search by type for the earlier chosen nuclide
                                search_string_type = '|'.join(['F', 'M', 'S'])
                            else:
                                # search by type for the earlier chosen nuclide
                                search_string_type = '|'.join([str(type_rad_elem)])

                            df_final = df[df['Type'].str.contains(search_string_type, na=False)]

                            if df_final.empty:
                                data_available = False
                                print("data for the specified type {} for radionuclide {} is not available. ". \
                                      format(type_rad_elem, rad))
                                if not df['Type'].empty:
                                    print("For {}, only following data is available at the moment: {}".format(rad, df))
                            else:
                                data_available = True
                                self.input_dict['type_rad'] = type_rad


                except Exception as ve:
                    print("input is not recognized")

        return type_rad, self.input_dict

    # option for considering progeny contribution in dose calculation: True or False
    def get_consider_progeny(self):
        consider_progeny = None
        while consider_progeny not in ['Y', 'N']:
            consider_progeny = input(
                color + "Do you want to consider the contribution of progeny for dose calculation? press Y/N \n")
            if consider_progeny == 'Y':
                self.input_dict['consider_progeny'] = True
                ignore_half_life = None
                while type(ignore_half_life) != int:
                    try:
                        ignore_half_life = ast.literal_eval(input(
                            color + "Provide half life value (in seconds) below which the effect of progenies will be ignored (Default: 1800). \n"))
                        if type(ignore_half_life) == int:
                            print(
                                "The half life value (in seconds) below which the effect of progenies will be ignored "
                                "is {} seconds".format(ignore_half_life))
                            self.input_dict['ignore_half_life'] = ignore_half_life
                        else:
                            print("Input must be an integer")

                    except Exception as ve:
                        print("Input must be an integer")

            elif consider_progeny == str('N'):
                self.input_dict['consider_progeny'] = False
            else:
                print("Input must be Y or N")
        return self.input_dict

    # release height of plume
    def get_release_height(self):
        release_height = None
        # while type(release_height) is None and not release_height.isnumeric():
        a = [False]
        while a:
            a = [_ for _ in [isinstance(release_height, int), isinstance(release_height, float)] if _ == True]
            try:
                release_height = ast.literal_eval(
                    input(color + "Provide release height of the plume (in metre) (e.g. 10): \n"))
                if isinstance(release_height, int) or isinstance(release_height, float):
                    print("release_height of the plume {} metre.".format(release_height))
                    self.input_dict['release_height'] = release_height
                    if release_height < 10:
                        print("N.B: No height-based correction factor will be applied for measurement height = 10 m "
                              "and release height < 10 m. With release height < 10 m and measurement height (M.H) >= "
                              "10 m, the corrected speed will be corrected considering measurement height of 10 m.")
                else:
                    a = [False]
                    print("input for release height must be a NUMBER.")
            except Exception as ve:
                a = [False]
                print("input for release height must be an integer.")

        return release_height, self.input_dict

    # list of downwind distance for calculation; must be a list e.g. [100,200] in metre[100,200,300,500,800,1000,1600,2000,3000,5000]
    def get_downwind_distances(self):
        downwind_distances = None
        while downwind_distances == None or type(downwind_distances) != list or \
                len([each for each in downwind_distances if type(each) == int]) != len(downwind_distances):
            try:
                downwind_distances = ast.literal_eval(
                    input(color + "Provide the list of downwind distances (in metre) for "
                                  "dilution factor/dose calculation (e.g. [100,200] \n"))
                if type(downwind_distances) == list and \
                        len([each for each in downwind_distances if type(each) == int]) == len(downwind_distances):
                    print(
                        "downwind distances (in metre) for which dose to be calculated: {}".format(downwind_distances))
                    self.input_dict['downwind_distances'] = downwind_distances
                else:
                    print('input not recognized. desired format: [100,200]')
            except Exception as ve:
                print('input not recognized. desired format: [100,200]')

        return downwind_distances, self.input_dict

    # the total dose (component wise i.e. total inhalation, total ground shine, total submersion and total ingestion dose)
    # will be printed for the plant boundary (metre); note that any distance from the input list
    # of "downwind_distances" can be designated plant_boundary
    def get_plant_boundary(self, downwind_distances):
        plant_boundary = None
        while type(plant_boundary) != int or plant_boundary not in downwind_distances:
            try:
                plant_boundary = ast.literal_eval(input(
                    color + "Provide the distance to the plant boundary (metre) for printing the contribution of total dose "
                            "(inhalation, ground shine, submersion and ingestion dose) at the plant boundary (metre). "
                            "\nNote that user must choose a distance from the list of downwind distances provided earlier. \n"))

                if type(plant_boundary) == int:
                    print("plant boundary: {} metre".format(plant_boundary))
                    self.input_dict['plant_boundary'] = plant_boundary
                else:
                    print("input for plant_boundary must be an integer.")
            except Exception as ve:
                print(
                    "input for plant_boundary must be an integer and must be from the provided list of downwind distances.")
        return plant_boundary, self.input_dict

    def get_annual_discharge_bq_rad_list(self, rads_list):
        annual_discharge_bq_rad_list = None
        while annual_discharge_bq_rad_list is None or not isinstance(annual_discharge_bq_rad_list, list) \
                or len(annual_discharge_bq_rad_list) != len(rads_list) \
                or any(not isinstance(each, (int, float)) for each in annual_discharge_bq_rad_list):
            try:
                annual_discharge_bq_rad_list = ast.literal_eval(
                    input("Provide the annual discharges for each selected radionuclides in Bq/year "
                          "(e.g. [1000000, 1000000, 2000000]): \n"))
                if isinstance(annual_discharge_bq_rad_list, list) \
                        and len(annual_discharge_bq_rad_list) == len(rads_list) \
                        and all(isinstance(each, (int, float)) for each in annual_discharge_bq_rad_list):
                    print("The annual discharge values {}".format(annual_discharge_bq_rad_list))
                    self.input_dict['annual_discharge_bq_rad_list'] = annual_discharge_bq_rad_list
                else:
                    annual_discharge_bq_rad_list = None
                    print("Please provide the annual discharge values in proper format "
                          "\n(e.g. [1000000, 1000000, 2000000] if three radionuclides are chosen)")
            except Exception as ve:
                print("Please provide the annual discharge values in proper format "
                      "\n(e.g. [1000000, 1000000, 2000000] if three radionuclides are chosen)")
        return annual_discharge_bq_rad_list, self.input_dict

    # time of daily operation:
    # follows 24-hour time format. 0-24 indicate 24 hour operation. This is only used for computation of dilution factor. Based on the
    # user specific data, only relevant met data for the mentioned period is used to compute dilution factor.
    # this is useful for the cases where the plant is operating during specific time of the day.
    def get_operation_time(self):
        operation_time = None
        while (type(operation_time) != list) \
                or ([0 <= operation_time[0] <= 24 and 0 <= operation_time[1] <= 24][0] != True) \
                or (type(operation_time) != list) \
                or len([each for each in operation_time if type(each) == int]) != len(operation_time):
            try:
                operation_time = ast.literal_eval(
                    input(color + "Provide the operation start and end time in a day (e.g. "
                                  "[0, 24] ;it indicates plant runs 24 hours; input [8,"
                                  "20] indicates plant start at 8 AM and shuts down at 8 "
                                  "PM) \n"))
                if type(operation_time) == list and len([each for each in operation_time if type(each) == int]) == len(
                        operation_time):
                    start_operation_time = operation_time[0]
                    end_operation_time = operation_time[1]
                    self.input_dict['start_operation_time'] = start_operation_time
                    self.input_dict['end_operation_time'] = end_operation_time
                else:
                    print("Please provide the operation start and end time in a day in proper format (e.g. [8,10])."
                          "\nThe value should be in the range 0 to 24")

            except Exception as ve:
                print(
                    "Please provide the operation start and end time in a day in proper format (e.g. [8,10]). The value "
                    "should be in the range 0 to 24")

        return operation_time, self.input_dict

    # measurement height of meteorological data
    def get_measurement_height(self):
        measurement_height = None
        while type(measurement_height) != int:
            measurement_height = input(color + "Provide measurement_height at which metereological data is measured. "
                                               "In case metereological data not available, a hypothetical measurement "
                                               "height should be given referring to the wind speed for conservative computation "
                                               "N.B: Wind speeds measured at the 10 m level should be used for those "
                                               "times when the effluent plume is considered to be a ground level "
                                               "release.)  (input example: 10): \n")

            try:
                if type(ast.literal_eval(measurement_height)) == int and ast.literal_eval(measurement_height) >= 10:
                    measurement_height = ast.literal_eval(measurement_height)
                    print(
                        "measurement height at which metereological data is measured is {} metre".format(
                            measurement_height))
                    self.input_dict['measurement_height'] = measurement_height

                else:
                    print("input for measurement_height must be an integer and the value should be equal or more than "
                          "10 m.")

            except Exception as ve:
                print("input for measurement_height must be an integer.")

        return measurement_height, self.input_dict

    def get_calm_correction(self):
        calm_correction = None
        while calm_correction not in ['Y', 'N']:
            calm_correction = input("Do you want to incorporate calm correction? press Y or N \n")
            if calm_correction == 'Y':
                print("calm correction will be incorporated")
                self.input_dict['calm_correction'] = True
            elif calm_correction == 'N':
                print("calm correction won't be incorporated")
                self.input_dict['calm_correction'] = False
            else:
                print('input is not recognized.')
        return calm_correction, self.input_dict

    def get_run_dose_computation(self):
        run_dose_computation = None
        while run_dose_computation not in ['Y', 'N']:
            run_dose_computation = input("Do you want to perform DOSE computation? press Y or N \n")
            if run_dose_computation == 'Y':
                print("DOSE Computation will be performed.")
                self.input_dict['run_dose_computation'] = True
            elif run_dose_computation == 'N':
                print("DOSE Computation won't be performed.")
                self.input_dict['run_dose_computation'] = False
            else:
                print('input is not recognized.')
        return run_dose_computation, self.input_dict

    # provide consumption time of food for instantaneous release; default is 30 days
    # unit is in days
    def get_consumption_time_food(self):
        consumption_time_food = None
        while type(consumption_time_food) != int:
            try:
                consumption_time_food = input(color + "Provide the consumption time (unit days) of food "
                                                      "(default is 30 days) \n(e.g. input: 30) \n")
                if not consumption_time_food:
                    consumption_time_food = 30
                    self.input_dict['consumption_time_food'] = consumption_time_food
                    print('For ingestion dose: time frame to consume food: {} days'.format(consumption_time_food))
                elif type(ast.literal_eval(consumption_time_food)) == int:
                    consumption_time_food = ast.literal_eval(consumption_time_food)
                    self.input_dict['consumption_time_food'] = consumption_time_food
                    print('For ingestion dose: time frame to consume food: {} days'.format(consumption_time_food))
                else:
                    print("input for the time frame must be an integer.")
            except Exception as ve:
                print("input for the time frame must be an integer.")
        return get_consumption_time_food, self.input_dict

    def get_instantaneous_release_bq_list(self, rads_list):
        instantaneous_release_bq_list = None
        while instantaneous_release_bq_list is None or not isinstance(instantaneous_release_bq_list, list) \
                or len(instantaneous_release_bq_list) != len(rads_list) \
                or any(not isinstance(each, (int, float)) for each in instantaneous_release_bq_list):
            try:
                instantaneous_release_bq_list = ast.literal_eval(input(
                    "Provide the instantaneous/accidental discharges for each selected radionuclides "
                    "in Bq/year (e.g. [1000000, 1000000, 2000000]): \n"))

                if isinstance(instantaneous_release_bq_list, list) \
                        and len(instantaneous_release_bq_list) == len(rads_list) \
                        and all(isinstance(each, (int, float)) for each in instantaneous_release_bq_list):
                    print("The instantaneous/accidental discharge values: {}".format(instantaneous_release_bq_list))
                    self.input_dict['instantaneous_release_bq_list'] = instantaneous_release_bq_list
                else:
                    instantaneous_release_bq_list = None
                    print('Input not recognized. Desired format: [100, 200]')
            except Exception as ve:
                print('Input not recognized. Desired format: [100, 200]')
        return instantaneous_release_bq_list, self.input_dict

    def get_max_dilution_factor(self, downwind_distances):
        max_dilution_factor = None
        while type(max_dilution_factor) != dict:
            try:
                max_dilution_factor = input(color + "Provide the value of dilution factor in dict format"
                                                    "\n(e.g. {dist_1: max_dilution_factor1, dist_2: "
                                                    "max_dilution_factor2}) \n")

                if type(ast.literal_eval(max_dilution_factor)) == dict \
                        and len(ast.literal_eval(max_dilution_factor)) == len(downwind_distances):
                    max_dilution_factor = ast.literal_eval(max_dilution_factor)
                    print(
                        "User provide dilution factor {} will be used for dose computation".format(max_dilution_factor))
                    self.input_dict['list_max_dilution_factor'] = max_dilution_factor

                else:
                    print("input for dilution factor must be a dictionary.")
            except Exception as ve:
                print("input for dilution factor must be a dictionary."
                      "The maximum dilution factor across the stability category for each specified downwind distances "
                      "must be given in dict format.")
        return max_dilution_factor, self.input_dict

    def get_weathering_corr_gs(self):
        weathering_corr = None
        while weathering_corr not in ['Y', 'N']:
            weathering_corr = input(
                color + "Do you like to incorporate weathering correction factor for ground-shine dose computation? "
                        "Press Y or N \n")
            if weathering_corr == 'Y':
                self.input_dict['weathering_corr'] = True
            elif weathering_corr == 'N':
                self.input_dict['weathering_corr'] = False
            else:
                print("The input should be either Y or N")
        return self.input_dict

    def get_weathering_corr_ingestion(self):
        weathering_corr = None
        while weathering_corr not in ['Y', 'N']:
            weathering_corr = input(
                color + "Do you like to incorporate weathering correction factor for ingestion dose computation? For "
                        "more details on this correction, refer to section 5.1.1.2 and Table VIII (page 63-64) IAEA "
                        "SRS 19\n. Press Y or N \n")
            if weathering_corr == 'Y':
                self.input_dict['weathering_corr_ingestion'] = True
            elif weathering_corr == 'N':
                self.input_dict['weathering_corr_ingestion'] = False
            else:
                print("The input should be either Y or N")
        return self.input_dict

    def get_n_value_plume_shine_dose(self):
        n = None
        while n not in [1, 2, 3]:
            n = ast.literal_eval(
                input(color + "provide the value of n for finding limit for triple integral evaluation. "
                              "The value should be within 1 and 3. "
                              "Higher value of n slightly increase the accuracy but significantly increase the computaional cost "
                              "(Default value is n=1) \n"))
            if 1 <= n <= 3:
                # for plume shine dose calculation. default=1, For high accuracy use n = 2 or 3, it exponentially increases the computational time
                print("chosen value of n for plume shine dose computation is {}".format(n))
                self.input_dict['n'] = int(n)
            else:
                print('the input must be an integer and between 1 to 3')
        return n, self.input_dict

    def get_inges_param_dict_adult(self):
        inges_param_dict_option_adult = None
        inges_param_dict_adult = {'DID_veg': 1.050, 'DID_milk': 0.500, 'DID_meat': 0.040, 'DID_fish': 0.050,
                                  'DID_water_and_beverage': 0.002}
        while (inges_param_dict_option_adult != 'D' and type(inges_param_dict_option_adult) != dict):
            try:
                inges_param_dict_option_adult = input(color + "Regarding ingestion dose computation for ADULT, "
                                                              "default values of the following parameters are printed value. "
                                                              "\n# additional parameters used for ingestion dose calculation (Adult); if not provided default values will be used and will be printed in log file. "
                                                              "\n# DID_veg = annual intake of vegetables # kg/day "
                                                              "\n# DID_milk = annual intake of milk # l/day "
                                                              "\n# DID_meat = annual intake of meat # kg/day "
                                                              "\n# DID_fish = annual intake of fish # kg/day "
                                                              "\n# DID_water_and_beverage = annual intake of water_and_beverages # m3/day "
                                                              "\n"
                                                              "\ninges_param_dict_adult: {DID_veg: 1.050, DID_milk: 0.500, DID_meat: 0.040, "
                                                              "DID_fish: 0.050, DID_water_and_beverage: 0.002} "
                                                              "In case you like to change value of any parameter, "
                                                              "please provide the name of the parameter with its desired "
                                                              "value in following format: {'DID_veg': 3, 'DID_milk': 4}. "
                                                              "\n If you like to use default value press D. \n")

                if inges_param_dict_option_adult == 'D':
                    print(
                        "for ADULT DOSE COMPUTATION: following values will be used for ingestion dose calculation: {}".format(
                            inges_param_dict_adult))
                    self.input_dict['inges_param_dict_adult'] = inges_param_dict_adult
                    inges_param_dict_option_adult = ast.literal_eval(str(inges_param_dict_option_adult))

                elif type(ast.literal_eval(str(inges_param_dict_option_adult))) == dict:
                    inges_param_dict_option_adult = ast.literal_eval(str(inges_param_dict_option_adult))
                    for k, v in inges_param_dict_option_adult.items():
                        inges_param_dict_adult[k] = v
                    self.input_dict['inges_param_dict_adult'] = inges_param_dict_adult
                    inges_param_dict_option_adult = ast.literal_eval(str(inges_param_dict_option_adult))
                    print(
                        "for ADULT DOSE COMPUTATION: following values will be used for ingestion dose calculation: {}".format(
                            inges_param_dict_adult))
                else:
                    print("please follow the instructions for providing input.")

            except Exception as ve:
                print("please follow the instructions for providing input for adult dose computation.")
        return inges_param_dict_adult, self.input_dict

    def get_inges_param_dict_infant(self):
        inges_param_dict_infant = {'DID_veg': 0.215, 'DID_milk': 0.400, 'DID_meat': 0.0032876, 'DID_fish': 0.004109,
                                   'DID_water_and_beverage': 0.0007123}

        inges_param_dict_option_infant = None
        while (inges_param_dict_option_infant != 'D' and type(inges_param_dict_option_infant) != dict):
            try:
                inges_param_dict_option_infant = input(color + "Regarding ingestion dose computation for INFANT, "
                                                               "default values of the following parameters are printed value. "
                                                               "\n# additional parameters used for ingestion dose calculation (infant); if not provided default values will be used and will be printed in log file "
                                                               "\n# DID_veg = annual intake of vegetables # kg/day "
                                                               "\n# DID_milk = annual intake of milk # l/day "
                                                               "\n# DID_meat = annual intake of meat # kg/day "
                                                               "\n# DID_fish = annual intake of fish # kg/day "
                                                               "\n# DID_water_and_beverage = annual intake of water_and_beverages # m3/day "
                                                               "\n"
                                                               "\ninges_param_dict_infant: {DID_veg: 0.215, DID_milk: 0.400, DID_meat: 0.0032876, DID_fish: 0.004109, DID_water_and_beverage: 0.0007123} "
                                                               "\nIn case you like to change value of any parameter, "
                                                               "please provide the name of the parameter with its desired "
                                                               "value in following format: {'DID_veg': 3, 'DID_milk': 4} "
                                                               "\n If you like to use default value press D. ")

                if inges_param_dict_option_infant == 'D':
                    print(
                        "for INFANT DOSE COMPUTATION: following values will be used for ingestion dose calculation: {}".format(
                            inges_param_dict_infant))
                    self.input_dict['inges_param_dict_infant'] = inges_param_dict_infant
                    inges_param_dict_option_infant = ast.literal_eval(str(inges_param_dict_option_infant))
                    break

                elif type(ast.literal_eval(str(inges_param_dict_option_infant))) == dict:
                    inges_param_dict_option_infant = ast.literal_eval(str(inges_param_dict_option_infant))
                    for k, v in inges_param_dict_option_infant.items():
                        inges_param_dict_infant[k] = v
                    self.input_dict['inges_param_dict_infant'] = inges_param_dict_infant
                    print(
                        "for INFANT DOSE COMPUTATION: following values will be used for ingestion dose calculation: {}".format(
                            inges_param_dict_infant))
                    break
                else:
                    print("please follow the instructions for providing input.")
            except Exception as ve:
                print("please follow the instructions for providing input for infant dose computation.")

        return inges_param_dict_infant, self.input_dict

    def get_veg_type_list_h3c14(self):
        veg_type_list = None
        while veg_type_list not in ['leafy_vegetables', 'non_leafy_vegetables', 'root_crops', 'all_others']:
            veg_type_list = input(color + "Provide the terrestrial vegetable option for ingestion dose computation "
                                          "of tritium/C-14. "
                                          "\nOptions: leafy_vegetables,non_leafy_vegetables, root_crops, all_others."
                                          "\n(e.g. input: leafy_vegetables) \n")
            if veg_type_list not in ['leafy_vegetables', 'non_leafy_vegetables', 'root_crops', 'all_others']:
                print("cannot identify the terrestrial vegetable option.")
            else:
                print("Terrestrial vegetable for tritium dose computation {}".format(veg_type_list))
                self.input_dict['veg_type_list'] = veg_type_list
        return veg_type_list, self.input_dict

    def get_climate_h3c14(self):
        # options: arctic, Mediterranean, meritime, Continental
        climate = None
        while climate not in ['arctic', 'Mediterranean', 'meritime', 'Continental']:

            climate = input(color + "Provide the climate condition for tritium dose computation"
                                    "\n(Options: arctic, Mediterranean, meritime, Continental)"
                                    "\n(e.g. Continental) \n")

            if climate not in ['arctic', 'Mediterranean', 'meritime', 'Continental']:
                print(
                    "please choose climate condition among these four options: arctic, Mediterranean, meritime, Continental")
            else:
                print("chosen climate is {}".format(climate))
                self.input_dict['climate'] = climate
        return climate, self.input_dict

    def get_animal_feed_type_h3c14(self):
        animal_feed_type = None
        while animal_feed_type not in ['leafy_vegetables', 'non_leafy_vegetables', 'root_crops', 'all_others']:
            animal_feed_type = input(color + "Provide the list of feed type to animals for ingestion dose computation "
                                             "of tritium/C-14. "
                                             "\nOptions: leafy_vegetables,non_leafy_vegetables, root_crops, all_others."
                                             "\n(e.g. input: leafy_vegetables) \n")

            if animal_feed_type not in ['leafy_vegetables', 'non_leafy_vegetables', 'root_crops', 'all_others']:
                print("can not identify the animal feed type ")
            else:
                print("animal feed type for tritium dose computation {}".format(animal_feed_type))
                self.input_dict['animal_feed_type'] = animal_feed_type
        return animal_feed_type, self.input_dict

    def get_animal_product_list_for_tritium(self):
        animal_product_list_for_tritium = None
        while type(animal_product_list_for_tritium) != list:
            try:
                animal_product_list_for_tritium = ast.literal_eval(
                    input(color + "For tritium: Provide the list of animal product for ingestion dose computation. "
                                  "\nOptions: cow_milk,goat_milk, goat_meat,lamb_meat,beef_meat,pork_meat, broiler_meat, egg. "
                                  "\n(e.g. input for cow milk and goat meat is ['cow_milk', 'goat_meat']) \n"))

                if type(animal_product_list_for_tritium) == list:
                    print(
                        "animal product list for tritium dose computation {}".format('animal_product_list_for_tritium'))
                    self.input_dict['animal_product_list_for_tritium'] = animal_product_list_for_tritium
                else:
                    print(
                        "can't recognize input (not a list) (input formar should follow this: ['cow_milk', 'goat_meat']).")

            except Exception as ve:
                print(
                    "can't recognize input (not a list) (input formar should follow this: ['cow_milk', 'goat_meat']).")

        return animal_product_list_for_tritium, self.input_dict

    def get_animal_product_list_for_C14(self):
        animal_product_list_for_C14 = None
        while type(animal_product_list_for_C14) != list:
            try:
                animal_product_list_for_C14 = ast.literal_eval(
                    input(color + "For C-14: Provide the list of animal product for ingestion dose computation. "
                                  "\nOptions: cow_milk,goat_milk, goat_meat,lamb_meat,beef_meat,pork_meat, broiler_meat, egg. "
                                  "\n(e.g. input for cow milk and goat meat is ['cow_milk', 'goat_meat']) \n"))

                if type(animal_product_list_for_C14) == list:
                    print("animal product list for C14 dose computation {}".format('animal_product_list_for_C14'))
                    self.input_dict['animal_product_list_for_C14'] = animal_product_list_for_C14
                else:
                    print(
                        "can't recognize input (not a list) (input formar should follow this: ['cow_milk', 'goat_meat']).")
            except Exception as ve:
                print(
                    "can't recognize input (not a list) (input formar should follow this: ['cow_milk', 'goat_meat']).")

        return animal_product_list_for_C14, self.input_dict

    def get_soiltype(self):
        soiltype = None
        # for ingestion dose computation
        while soiltype not in ['peatsoil', 'othersoil']:
            try:
                soiltype = ast.literal_eval(
                    input(color + 'Provide the soil type for ingestion dose computation. Options are '
                                  'either "peatsoil" or "othersoil" \n'))
                if soiltype in ['peatsoil', 'othersoil']:
                    print("chosen soil type is {}".format(soiltype))
                    self.input_dict['soiltype'] = str(soiltype)
                else:
                    print("suggested soil type not supported. Ensure the soil name is provided within quotation.")
            except Exception as ve:
                print("Input not recognized")
        return soiltype, self.input_dict

    def get_inges_param_dict(self):
        inges_param_dict_option = None
        inges_param_dict = {'alpha_wet_crops': 0.3, 'alpha_dry_forage': 3, 't_e_food_crops': 60, 't_e_forage_grass': 30,
                            't_b': 11000, 't_h_wet_crops': 14, 't_h_animal_pasture': 0, 't_h_animal_stored_feed': 90,
                            'C_wi': 0, 'f_p': 0.7, 'alpha': 3, 't_e': 30, 't_m': 1, 't_f': 20, 'q_m': 16, 'q_w': 0.06,
                            'q_f': 1.2, 'q_w_meat': 0.004}

        while inges_param_dict_option != 'D' and type(inges_param_dict_option) != dict:
            try:
                inges_param_dict_option = input(color + "For ingestion dose computation, "
                                                        "default values of the following parameters are printed value. "
                                                        "\n# alpha (m2/kg) = fraction of deposited activity intercepted by the edible portion of vegetation per unit mass; source: table VIII, page 64 of SRS 19 "
                                                        "\n# t_e (days) = the period that the crops are exposed to contamination during growing season. "
                                                        "\n# F_v (Bq/kg dry soil) = concentration for uptake of radionuclide from soil by edible parts of crops. "
                                                        "\n# t_b = duration of discharge of material in a day; for 30 years it is 11000. "
                                                        "\n# t_h_wet_crops (days) = delay time i.e. time interval between harvest and consumption of food. "
                                                        "\n# t_h_animal_pasture, t_h_animal_stored_feed = time in animal feed. "
                                                        "\n# f_p = fraction of the year that animal consumes fresh pasture vegetation. "
                                                        "\n# c_pi = concentration of radionuclide in stored feeds (Bq/Kg). "
                                                        "\n# t_h_delay_time = 90 (day). "
                                                        "\n# t_m (day) = average time between collection and human consumption of milk. "
                                                        "\n# t_f (day) = average time between collection and human consumption of meat. "
                                                        "\n# q_m (kg/day) = amount of feed (in dry matter) consumed per day. "
                                                        "\n# q_w (m3/d) = amount of water consumed by animal per day. Source: Table B1 B2 page 66 (ECPDA) "
                                                        "\n# C_wi (Bq/m3) = concentration of radionuclide in water."
                                                        "\n# q_f (Kg/d) = amount of feed consumed by animal; goat and sheep; meat producing animal."
                                                        "\n# q_w_meat (m3/day) = water intake of meat producing animal. and the associated default values of these parameters are: "
                                                        "\n"
                                                        "\ninges_param_dict = {alpha_wet_crops: 0.3, alpha_dry_forage: 3, t_e_food_crops: 60, "
                                                        "t_e_forage_grass: 30, t_b: 11000, t_h_wet_crops: 14, t_h_animal_pasture: 0, "
                                                        "t_h_animal_stored_feed: 90, C_wi: 0, f_p: 0.7, alpha: 3, t_e: 30, t_m: 1, "
                                                        "t_f: 20, q_m: 16, q_w: 0.06, q_f: 1.2, q_w_meat: 0.004}"
                                                        "In case you like to change value of any parameter, please provide the name of the parameter with its desired "
                                                        "value in following format: {'alpha_wet_crops': 0.1, 't_f': 4}. If you like to use default value press D \n")

                if inges_param_dict_option == 'D':
                    print("following values will be used for ingestion dose calculation: {}".format(inges_param_dict))
                    self.input_dict['inges_param_dict'] = inges_param_dict
                    inges_param_dict = ast.literal_eval(str(inges_param_dict_option))

                elif type(ast.literal_eval(str(inges_param_dict_option))) == dict:
                    inges_param_dict_option = ast.literal_eval(str(inges_param_dict_option))
                    for k, v in inges_param_dict_option.items():
                        inges_param_dict[k] = V
                    self.input_dict['inges_param_dict'] = inges_param_dict
                    print("following values will be used for ingestion dose calculation: {}".format(inges_param_dict))
                else:
                    print("please follow the instructions for providing input.")
            except Exception as ve:
                print("please follow the instructions for providing input.")
        return inges_param_dict, self.input_dict

    def get_run_plume_shine_dose(self):
        # for plume shine dose computation
        run_plume_shine_dose = None
        while run_plume_shine_dose not in ['Y', 'N']:
            run_plume_shine_dose = input(color + "Do you like to perform Plume shine dose computation? Press Y or N \n")
            if run_plume_shine_dose == 'Y':
                self.input_dict['run_plume_shine_dose'] = True
            elif run_plume_shine_dose == 'N':
                self.input_dict['run_plume_shine_dose'] = False
        return self.input_dict

    def get_num_days(self, excel_sheet_name):

        num_days = None
        # number of days in each year of met data
        while num_days == None or type(num_days) != list or \
                len([each for each in num_days if type(each) == int]) != len(num_days) \
                or len(num_days) != len(excel_sheet_name):
            try:
                num_days = ast.literal_eval(input(
                    color + "Provide the number of days of each year for which the metereological data is available (e.g. [365,365,365,366,365]): \n"))
                if type(num_days) == list and len([each for each in num_days if type(each) == int]) == len(num_days):
                    print('number of days for which metereological data available: {}'.format(num_days))
                    self.input_dict['num_days'] = num_days
                    # plot dilution factor sector wise, individual plots will be generated for specific down wind distances.
                    self.input_dict['plot_dilution_factor'] = True
                    # Calm correction
                    # sampling time for meteorological data; default= 60 minute; sampling_time correction is on TODO list
                    sampling_time = 60
                    self.input_dict['sampling_time'] = 60


                else:
                    print("Please provide the number of days of each year "
                          "for which the metereological data is available in proper format (e.g. [365,365,365,366,365]). "
                          "Ensure that the number of excel sheets should be equal to number of entries in the list. "
                          "\nFor example: [365, 366] indicates that two excel sheet present in which "
                          "first one has 365 day data, while second one has 366 days data")
            except Exception as ve:
                print("Please provide the number of days of each year "
                      "for which the metereological data is available in proper format (e.g. [365,365,365,366,365]). "
                      "Ensure that the number of excel sheets should be equal to number of entries in the list. "
                      "\nFor example: [365, 366] indicates that two excel sheet present in which "
                      "first one has 365 day data, while second one has 366 days data")
        return num_days, self.input_dict
        # num_days, input_dict = InpGenFunc.get_num_days(excel_sheet_name)
        # measurement height of meteorological data
        # measurement_height, input_dict = InpGenFunc.get_measurement_height()
        # calm correction (Yes/No)
        # calm_correction, input_dict = InpGenFunc.get_calm_correction()

    def get_column_names(self):
        column_names = None
        # column names of excel sheet containing meteorological data. The sequence/order of columns should be strictly followed.
        # for instance, hour should be the first index followed by speed, direction and stability category
        while type(column_names) != list:
            try:
                column_names = ast.literal_eval(input(
                    color + "provide the name of four columns (time, wind speed, direction and stability class) in excel sheets for accessing metereological data (e.g. ['HOUR', 'WS 10m(kmph)', 'DIR at 10m', 'STBCLASS']): \n"))
                if type(column_names) == list:
                    print("name of the columns within sheet of excel file: {}".format(column_names))
                    self.input_dict['column_names'] = column_names

                else:
                    print(
                        "Please provide the name of four columns  in proper format (e.g.  ['HOUR', 'WS 10m(kmph)', 'DIR at 10m', 'STBCLASS'])")
            except Exception as ve:
                print(
                    "Please provide the name of four columns  in proper format (e.g.  ['HOUR', 'WS 10m(kmph)', 'DIR at 10m', 'STBCLASS'])")

        return column_names, self.input_dict

    def get_excel_sheet_name(self):
        excel_sheet_name = None
        # sheet name; it is anticipated that year of data is the sheet name; multiple year may be provided; ['2017','2018']
        while type(excel_sheet_name) != list:
            try:
                excel_sheet_name = ast.literal_eval(input(color + "provide the name of excel sheets for accessing meteorological data \
                                            (e.g. ['2017','2018'] or ['sheet1', 'sheet2']): \n"))
                if type(excel_sheet_name) == list:
                    print("name of the excel sheets are: {}".format(excel_sheet_name))
                    self.input_dict['excel_sheet_name'] = excel_sheet_name


            except Exception as ve:
                print(
                    "Please provide the name of sheets  in proper format (e.g. ['2017','2018'] or ['sheet1', 'sheet2'])")

        return excel_sheet_name, self.input_dict

    def get_details_have_dilution_factor(self):
        have_dilution_factor = None
        while have_dilution_factor not in ['Y', 'N']:
            try:
                have_dilution_factor = input(
                    color + "Do you like to provide precomputed dilution factor for the dose computation. Press Y or N \n")
                if have_dilution_factor == 'N':
                    self.input_dict['max_dilution_factor'] = None
                    self.input_dict['have_dilution_factor'] = False

                elif have_dilution_factor == 'Y':
                    self.input_dict['have_dilution_factor'] = True

                else:
                    print("input must be Y or N.")
            except Exception as ve:
                print("input must be Y or N.")

        return have_dilution_factor, self.input_dict

    def get_release_scenario(self):
        # Type of release scenario
        release_scenario = None
        while release_scenario not in [1, 2]:
            try:
                release_scenario = int(
                    input(color + "press 1 for Long Term Continuous Release, 2 for Instantenous Release: \n"))
                if release_scenario == 1:
                    print("you opted for Long Term Continuous Release")
                    self.input_dict['long_term_release'] = True
                    self.input_dict['single_plume'] = False

                # RELEASE SCENARIO 2; INSTANTENOUS RELEASE
                elif release_scenario == 2:
                    print("you opted for Instantenous Release ")
                    self.input_dict['long_term_release'] = False
                    self.input_dict['single_plume'] = True

                else:
                    # raise Exception('Please press either 1 or 2.')
                    print("Please press either 1 or 2.")

            except Exception as ve:
                print("Please press either 1 or 2.")
        return release_scenario, self.input_dict

    def get_master_eq_single_plume(self):
        # for plume shine dose computation
        max_conc_plume_central_line_gl = None
        while max_conc_plume_central_line_gl not in ['Y', 'N']:
            max_conc_plume_central_line_gl = input(color + "Do you like to perform dose computation that"
                                                           " projects maximum concentration along"
                                                           " the projection of plume centerline on the ground (i.e. Y = 0, Z=0)? Press Y or N \n")
            if max_conc_plume_central_line_gl == 'Y':
                self.input_dict['max_conc_plume_central_line_gl'] = True
                self.input_dict['Y'] = 0
                self.input_dict['Z'] = 0
            elif max_conc_plume_central_line_gl == 'N':
                self.input_dict['max_conc_plume_central_line_gl'] = False
                val_YZ, _ = self.get_YZ()
        return self.input_dict

    def get_master_eq_continuous_plume(self):
        # for plume shine dose computation
        max_conc_plume_central_line_gl = None
        while max_conc_plume_central_line_gl not in ['Y', 'N']:
            max_conc_plume_central_line_gl = input(color + "Do you like to compute ground level time-integrated "
                                                           "concentration (i.e. Z=0)? Press Y or N \n")
            if max_conc_plume_central_line_gl == 'Y':
                self.input_dict['max_conc_plume_central_line_gl'] = True
                self.input_dict['Z'] = 0
            elif max_conc_plume_central_line_gl == 'N':
                self.input_dict['max_conc_plume_central_line_gl'] = False
                val_Z, _ = self.get_Z()
        return self.input_dict

    def get_YZ(self):
        """
        ask Y and Z values
        @return: Y and Z values to be used in master equation
        """
        YZ = None
        while YZ is None or type(YZ) != list or len(YZ) != 2 or \
                len([each for each in YZ if type(each) == int]) != len(YZ):
            try:
                YZ = ast.literal_eval(input(
                    color + "Please provide the value of Y and Z in proper format (e.g. [0,10] for Y = 0 "
                            "and Z = 10 metre): \n"))
                if type(YZ) == list and len([each for each in YZ if type(each) == int]) == len(YZ) and len(YZ) == 2:
                    print('Dilution factor will be computed using Y = {} m and Z = {} m'.format(YZ[0], YZ[1]))
                    self.input_dict['Y'] = YZ[0]
                    self.input_dict['Z'] = YZ[1]
                else:
                    print(
                        "Please provide the value of Y and Z in proper format (e.g. [0,10] for Y = 0 and Z = 10) "
                        "metre). ")
            except Exception as ve:
                print("Please provide the value of Y and Z in proper format (e.g. [0,10] for Y = 0 and Z = 10 metre).")

        return YZ, self.input_dict

    def get_Z(self):
        """
        ask Y and Z values
        @return: Y and Z values to be used in master equation
        """
        Z = None
        while Z is None or type(Z) != list or len(Z) != 1 or \
                len([each for each in Z if type(each) == int]) != len(Z):
            try:
                Z = ast.literal_eval(input(
                    color + "Please provide the value of Z in proper format (e.g. [10] for "
                            " Z = 10 metre): \n"))
                if type(Z) == list and len([each for each in Z if type(each) == int]) == len(Z) and len(Z) == 1:
                    print('Dilution factor will be computed using Z = {} m'.format(Z[0]))
                    self.input_dict['Z'] = Z[0]

                else:
                    print(
                        "Please provide the value of Z in proper format (e.g. [0,10] for Z = 10) "
                        "metre). ")
            except Exception as ve:
                print("Please provide the value of Z in proper format (e.g. [10] for Z = 10 metre).")

        return Z, self.input_dict

    def get_age_group(self):
        # age list for calculation; must be a list e.g. [1,18]; 1 for infant and 18 for adult.
        age_group = None
        while age_group is None or type(age_group) != list or len(
                [each for each in age_group if type(each) == int]) != len(
            age_group):
            try:
                age_group = ast.literal_eval(input(color + "Provide the list of ages for dose computation "
                                                           "(e.g. [1,18]; 1 for infant and 18 for adult): \n"))
                self.input_dict['age_group'] = age_group
            except Exception as ve:
                print('The input for age group not recognized. Desired format is [1,5,18].')
        return age_group, self.input_dict

    def get_path_met_file(self):
        path_met_file = '.'
        while not os.path.isfile(path_met_file):
            try:
                # path of excel file containing meteorological data.
                path_met_file = ast.literal_eval(input(color + "Provide path of excel file "
                                                               "containing meteorological data of the site (e.g. 'met_data/Met_data_5Yr.xlsx'): \n"))

                if os.path.isfile(path_met_file):
                    print('provided path for metereological excel file: {}'.format(path_met_file))
                    self.input_dict['path_met_file'] = path_met_file
            except Exception as ve:
                print('file not found. please check the path.')
        return path_met_file, self.input_dict

    def get_have_met_data_release_scenario_1(self, release_scenario):
        have_met_data = None
        self.input_dict['plot_dilution_factor'] = False
        if release_scenario == 1:
            # for plume shine dose computation
            while have_met_data not in ['Y', 'N']:
                have_met_data = input(color + "Do you have meteorological data? Press Y or N \n")
                if have_met_data == 'Y':
                    self.input_dict['have_met_data'] = True
                    self.input_dict['scaling_dilution_factor_based_on_met_data_speed_distribution'] = True

                elif have_met_data == 'N':
                    self.input_dict['have_met_data'] = False
                    self.input_dict['max_dilution_factor'] = None
                    self.input_dict['scaling_dilution_factor_based_on_met_data_speed_distribution'] = False
                else:
                    print("input must be either 'Y' or 'N'")
        else:
            pass
        return self.input_dict['have_met_data'], self.input_dict

    def get_have_met_data(self, have_dilution_factor, release_scenario):
        if release_scenario == 1 or release_scenario == 2:
            if have_dilution_factor == 'Y':
                self.input_dict['have_met_data'] = False
            elif have_dilution_factor == 'N':
                self.input_dict['have_met_data'] = True
            else:
                print("The input must be either Y or N")
        else:
            pass
        return self.input_dict['have_met_data'], self.input_dict

    def get_have_met_data_release_scenario_2(self, release_scenario):
        have_met_data = None
        if release_scenario == 2:
            # for plume shine dose computation

            while have_met_data not in ['Y', 'N']:
                have_met_data = input(color + "Do you have meteorological data? Press Y or N \n")
                if have_met_data == 'Y':
                    self.input_dict['have_met_data'] = True
                    self.input_dict['scaling_dilution_factor_based_on_met_data_speed_distribution'] = True

                elif have_met_data == 'N':
                    self.input_dict['have_met_data'] = False
                    self.input_dict['max_dilution_factor'] = None
                    self.input_dict['scaling_dilution_factor_based_on_met_data_speed_distribution'] = False
                else:
                    print("input must be either 'Y' or 'N'")
        else:
            pass
        return have_met_data, self.input_dict

    def get_mean_speed_stab_cat_wise(self, release_scenario=2):
        ask_mean_speed_data = None
        if release_scenario == 2:
            while type(ask_mean_speed_data) != list:
                # for plume shine dose computation; single plume; for the scaling of the dilution factor
                try:
                    ask_mean_speed_data = ast.literal_eval(input(color + "Provide mean speed data for"
                                                                         "all six stability categories in array form. "
                                                                         "example: [1,1,2,4,5,1] \n"))

                    if type(ask_mean_speed_data) == list and len(ask_mean_speed_data) == 6:
                        print("provided mean speed data: {}".format(ask_mean_speed_data))
                        self.input_dict['ask_mean_speed_data'] = ask_mean_speed_data
                    else:
                        print("please provide an array of six mean speed associated "
                              "with the siz stability categories")

                except Exception as ve:
                    print("Please provide the mean speed of all stability categories in proper format"
                          " (e.g. [1,2, 4, 5, 1, 2]")
        else:
            pass

        return ask_mean_speed_data, self.input_dict

    def scale_dilution_factor_with_user_defined_mean_speed(self):
        # for single plume: ask user if they want to scale the dilution factor with user-define mean speed data.
        like_to_scale_with_mean_speed = None
        while like_to_scale_with_mean_speed not in ['Y', 'N']:
            like_to_scale_with_mean_speed = input(color + "Do you like to scale dilution factor with "
                                                          "your chosen mean speed? Press Y or N \n")
            if like_to_scale_with_mean_speed == 'Y':
                self.input_dict['like_to_scale_with_mean_speed'] = True

            elif like_to_scale_with_mean_speed == 'N':
                self.input_dict['like_to_scale_with_mean_speed'] = False
            else:
                print("input must be either 'Y' or 'N'")
        return like_to_scale_with_mean_speed, self.input_dict

    def get_exposure_period_for_ground_shine(self):
        exposure_period = None
        while not isinstance(exposure_period, (int, float)):
            try:
                exposure_period = input(color + "provide the exposure_period (year) for Ground shine dose computation. "
                                                "the generic value used in the SRS is 30 years.(refer to page 64, "
                                                "SRS 19, Table VIII). "
                                                "for decomissioning: 5 year may be choosen. \n ")
                if isinstance(ast.literal_eval(exposure_period), (int, float)):
                    # for plume shine dose calculation. default=1, For high accuracy use n = 2 or 3, it exponentially
                    # increases the computational time
                    exposure_period = ast.literal_eval(exposure_period)
                    print("chosen value of exposure_period for ground shine dose computation is {} years".format(
                        exposure_period))
                    self.input_dict['exposure_period'] = exposure_period
                else:
                    print('the input must be an number (unit: year)')
            except Exception as ve:
                print('Input is not recognized')

        return exposure_period, self.input_dict

    def get_output_file_name(self):
        # output file name
        output_file_name = None
        while output_file_name is None and type(output_file_name) != str:
            try:
                output_file_name = input(
                    color + "Provide the desired name of the output file (default: pydose_out): \n")
                if not output_file_name:
                    self.input_dict['output_file_name'] = 'pydose_out'
                    print("desired file name: pydose_out")
                else:
                    self.input_dict['output_file_name'] = output_file_name
                    print("desired file name: {}".format(output_file_name))
            except Exception as ve:
                print("desired file name format incorrect.")
        return output_file_name, self.input_dict

    def get_logdir_name(self):
        logdir_name = None
        while logdir_name is None and type(logdir_name) != str:
            logdir_name = input(color + "Provide the name of the directory for saving input, output and log file ("
                                        "default: loginfo): \n")
            if not logdir_name:
                self.input_dict['logdir_name'] = 'loginfo'
                print("Desired directory name: loginfo")
            else:
                self.input_dict['logdir_name'] = logdir_name
                print("Desired directory name: {}".format(logdir_name))
        return logdir_name, self.input_dict

    def get_input_file_name(self, logdir_name):
        input_file_name = None
        file_exists = None
        while type(input_file_name) != str:
            try:
                input_file_name = input(color + "Provide the desired name of the input file (default: auto_input): \n")
                # default
                if input_file_name == 'D' or input_file_name == '':
                    input_file_name = 'auto_input'
                    self.input_dict['input_file_name'] = input_file_name
                    print("desired input file name: auto_input")
                    # check if the file exist already
                    # don't exist
                    if not os.path.isfile('%s/%s.yaml' % (logdir_name, input_file_name)):
                        file_exists = False
                        # run auto input generator
                        print("________________________________________________________________________")
                        print("generated content of input file:", self.input_dict)
                        with open('%s/%s.yaml' % (logdir_name, input_file_name), 'w') as infile:
                            yaml.safe_dump(self.input_dict, infile, default_flow_style=True)
                        print("That's all we needed. Now sit back and relax while we perform the computation :).")
                        # subprocess.call("python main.py --config_file %s --output_file_name %s --logdir %s"%('auto_input.yaml', input_dict['output_file_name'], input_dict['logdir_name']), shell=True)

                    # file exists
                    else:
                        file_exists = True
                        print("________________________________________________________________________")
                        print(
                            "Input file already exists in the present directory and will be used for the computation.")


                else:
                    self.input_dict['input_file_name'] = input_file_name
                    print("desired file name: {}.yaml".format(input_file_name))
                    # check if the file exist already
                    # don't exist
                    if not os.path.isfile('%s/%s.yaml' % (logdir_name, input_file_name)):
                        file_exists = False
                        # run auto input generator
                        print("________________________________________________________________________")
                        # print("generated content of input file:", self.input_dict)
                        with open('%s/%s.yaml' % (logdir_name, input_file_name), 'w') as infile:
                            yaml.safe_dump(self.input_dict, infile, default_flow_style=True)
                        # print("That's all we needed. Now sit back and relax while we perform the computation :).")
                        # subprocess.call("python main.py --config_file %s --output_file_name %s --logdir %s"%('auto_input.yaml', input_dict['output_file_name'], input_dict['logdir_name']), shell=True)

                    # file exists
                    else:
                        file_exists = True
                        print("________________________________________________________________________")
                        print(
                            "Input file already exists in the present directory and will be used for the computation.")
                        # subprocess.call("python main.py --config_file %s --output_file_name %s --logdir %s"%('auto_input.yaml', input_dict['output_file_name'], input_dict['logdir_name']), shell=True)

            except Exception as ve:
                print('error:', ve)
                print("desired file name format incorrect.")
        return input_file_name, file_exists, self.input_dict

    def check_existing_input_entries(self):
        filename = self.input_dict['input_file_name']
        logdir_name = self.input_dict['logdir_name']

        with open(f'{logdir_name}/{filename}.yaml', 'r') as infile:
            data = yaml.safe_load(infile)
            # Collect lengths in a list
            lengths = [(len(data['element_list']), len(data['type_rad']), len(data['rads_list']))]
            # Append based on the release type
            if data.get('long_term_release'):
                lengths.append(len(data['annual_discharge_bq_rad_list']))
            else:
                lengths.append(len(data['instantaneous_release_bq_list']))

            # Flatten the list of lengths
            flat_lengths = [item for sublist in lengths for item in
                            (sublist if isinstance(sublist, tuple) else (sublist,))]

            print("Checking the input File:")

            # Check if all lengths are equal
            if len(set(flat_lengths)) > 1:
                raise ValueError(
                    "Please check the input file. The following entries (element_list, type_rad, rads_list, "
                    "annual_discharge_bq_rad_list/instantaneous_release_bq_list) must have the same length."
                )
            else:
                print("PASSED initial checking!!")

    def get_pickle_it(self):
        try:
            pickle_it = self.input_dict['pickle_it']
        except KeyError:
            pickle_it = False
            self.input_dict['pickle_it'] = pickle_it
        return pickle_it

    def get_run_ps_dose_parallel(self):
        try:
            run_ps_dose_parallel = self.input_dict['run_ps_dose_parallel']
        except KeyError:
            run_ps_dose_parallel = True
            self.input_dict['run_ps_dose_parallel'] = run_ps_dose_parallel
        return run_ps_dose_parallel
