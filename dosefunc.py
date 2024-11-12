import numpy as np
import os
import pandas as pd
from scipy import integrate, special
import logging
import time
from joblib import Parallel, delayed

os.getcwd()
np.random.seed(seed=3007)


class DoseFunc:
    """
    Represents a class for dose computation.

    Args:
        - device: Description of device.
        - config: Configuration parameters.
        - logdir: Log directory (default is None).

    Attributes:
        - effective_lambda_of_rads: Effective lambda of radionuclides.
        - speed_dist_array: Speed distribution array.
        - mean_speed_stab_cat_wise: Mean speed stability category wise.
        - lambda_of_rads: Lambda of radionuclides.
        - dilution_factor_sectorwise: Dilution factor sector-wise.
        - dcf_list_submersion_corr: List of dose conversion factors with submersion correction.
        - dcfs_gs_corr: Dose conversion factors with ground shine correction.
        - max_dilution_factor: Maximum dilution factor.
        - ingestion_dose_vals: Ingestion dose values.
        - ingestion_dose_tritium_vals: Ingestion dose values for tritium.
        - dcfs_ingestion: Dose conversion factors for ingestion.
        - sub_dose_per_rad: Submersion dose per radionuclide.
        - dcf_list_submersion: List of dose conversion factors for submersion.
        - gs_dose_per_rad: Ground shine dose per radionuclide.
        - dcfs_gs: Dose conversion factors for ground shine.
        - inhalation_dcf: Inhalation dose conversion factor.
        - sigma_z: Sigma Z value.
        - AY: AY value.
        - sigma_y: Sigma Y value.
        - calm_correction_factors: Calm correction factors.
        - TJFD_ALL_MISSING_CORR: Triple joint frequency distribution with missing correction.
        - TJFD_ALL: Triple joint frequency distribution.
        - MET_DATA: Meteorological data.
        - filepath: Path to Excel file containing meteorological data.
        - sheet_names: Sheet names in Excel file containing meteorological data.
        - sampling_time: Sampling time.
        - colnames: Column names in meteorological data.
        - num_days: Number of days for meteorological data.
        - measurement_height: Measurement height of meteorological data.
        - release_height: Release height of plume.
        - downwind_distances: Downwind distances.
        - annual_discharge_bq_rad_list: Annual discharge of radionuclides list.
        - rads_list: List of radionuclides.
        - element_list: List of elements.
        - age_group: Age group for dose calculation.

    Raises:
        - ValueError: If required parameters are not provided or if conflicting configurations are set.
    """

    def __init__(self, device, config, log_file_name, logdir=None):
        super().__init__()
        self.effective_lambda_of_rads = None
        self._device = device
        self.config = config
        self.logdir = logdir
        self.speed_dist_array = None
        self.mean_speed_stab_cat_wise = None
        self.lambda_of_rads = None
        self.dilution_factor_sectorwise = None
        self.dcf_list_submersion_corr = None
        self.dcfs_gs_corr = None
        if self.config['have_dilution_factor']:
            dict_max_dilution_factor = config['list_max_dilution_factor']
        else:
            self.max_dilution_factor = None
        self.ingestion_dose_vals = None
        self.ingestion_dose_tritium_vals = None
        self.dcfs_ingestion = None
        self.sub_dose_per_rad = None
        self.dcf_list_submersion = None
        self.gs_dose_per_rad = None
        self.dcfs_gs = None
        self.inhalation_dcf = None
        self.sigma_z = None
        self.AY = None
        self.sigma_y = None
        self.calm_correction_factors = None
        self.TJFD_ALL_MISSING_CORR = None
        self.TJFD_ALL = None
        self.MET_DATA = None

        if config['have_met_data']:
            self.filepath = config['path_met_file']

        if self.config['long_term_release'] and self.config['single_plume']:
            raise ValueError("Choose the mode of release: Instantaneous or Continuous (long-term).")

        if not self.config['long_term_release'] and not self.config['single_plume']:
            raise ValueError("Choose the mode of release: Instantaneous or Continuous (long-term).")

        if config['have_met_data']:
            if not self.filepath:
                raise ValueError("Must provide the path of Excel file containing meteorological data.")

            self.sheet_names = config['excel_sheet_name']
            if not self.sheet_names:
                raise ValueError("Please provide the correct name of Excel sheet containing meteorological data.")

            self.sampling_time = config['sampling_time']

            self.colnames = config['column_names']

            self.num_days = config['num_days']
            if not self.num_days:
                raise ValueError("Please provide the number of days for which meteorological data is provided.")

            self.measurement_height = config['measurement_height']
            if not self.measurement_height:
                raise ValueError("Please provide the measurement height at which meteorological data is measured.")

        self.release_height = config['release_height']
        if not self.release_height:
            raise ValueError("Please provide the release height of plume.")

        self.downwind_distances = config['downwind_distances']
        if self.config['plant_boundary'] not in self.downwind_distances:
            self.downwind_distances.append(self.config['plant_boundary'])

        if not self.downwind_distances:
            raise ValueError("Please provide the list of downwind distances (in metre).")

        if self.config['run_dose_computation']:
            self.rads_list = config['rads_list']
            if not self.rads_list:
                raise ValueError("Must provide the list of radionuclides for dose estimation.")
            if 'H-3' in self.rads_list:
                if not self.config['animal_product_list_for_tritium']:
                    raise ValueError("animal_product_list_for_tritium not given")
                if not self.config['animal_product_list_for_tritium']:
                    raise ValueError("animal_feed_type not given")
                if not self.config['veg_type_list']:
                    raise ValueError("veg_type_list not given")
                if not self.config['climate']:
                    raise ValueError("climate not given")

            if config['long_term_release']:
                self.annual_discharge_bq_rad_list = config['annual_discharge_bq_rad_list']
                if not self.annual_discharge_bq_rad_list:
                    raise ValueError("Please provide annual discharge info (Bq/year) for each radionuclide.")

            self.element_list = config['element_list']
            if not self.element_list:
                raise ValueError("Must provide the list of elements for dose estimation.")

            self.age_group = config['age_group']
            if not self.age_group:
                raise ValueError("Please provide the list of age for which dose needs to be calculated.")

    def inhalation_dose(self, X1, age=1, max_dilutfac_for_distance_secperm3=None):
        """
            Calculate the inhalation dose for a given release scenario.

            Args:
                X1 (float): Distance from the source in meters.
                age (int, optional): Age of the recipient in years. Default is 1.
                max_dilutfac_for_distance_secperm3 (float, optional): Maximum dilution factor for the given distance in
                    seconds per cubic meter. If not provided, it will be calculated based on the release scenario.

            Returns:
                list: A list containing the inhalation dose for each radionuclide in the scenario.

            Raises:
                ValueError: If the age of the recipient is not a number.
        """
        if self.max_dilution_factor is None and self.config['have_met_data']:
            self.file_preprocessing()
            self.met_data_to_tjfd()
            self.missing_correction()
            self.calm_correction_factor_calc()
            self.dilution_per_sector(X1)
        if self.config['have_dilution_factor']:
            max_dilutfac_for_distance_secperm3 = self.dict_max_dilution_factor[X1]
        if not self.config['have_dilution_factor'] and max_dilutfac_for_distance_secperm3 is None:
            dilution_factor_sectorwise = self.dilution_per_sector(X1)
            # print('dilution_factor_sectorwise:', dilution_factor_sectorwise)
            # max_dilutfac_for_distance_secperm3 = max(dilution_factor_sectorwise)
            max_dilutfac_for_distance_secperm3 = np.array(dilution_factor_sectorwise).T[0].max()
            # max_dilutfac_for_distance_secperm3 = self.max_dilution_factor

        # m3/s
        age = age
        if age > 1:
            breathing_rate = 8400 / (365 * 3600 * 24)
        elif age <= 1:
            breathing_rate = 1400 / (365 * 3600 * 24)
        else:
            raise ValueError('The age of recipient must be a number.')

        # inhalation_dose = max_dilutfac_for_distance_secperm3 * annual_discharge_ci_rad_list
        dcfs = self.inhalation_dcf_list(master_file="library/RadioToxicityMaster.xls",
                                        sheet_name="Inhalation CED Sv per Bq Public", age=age)

        tot_inhalation_dose = 0
        inhal_dose_per_rad = []

        # for single_plume discharge is in Bq unit (Time integrated)
        if self.config['single_plume']:
            discharge_q = self.config['instantaneous_release_bq_list']
        # for continuous release, discharge in Bq/year unit
        else:
            discharge_q = self.annual_discharge_bq_rad_list

        for ndx, (rad, dcf) in enumerate(zip(discharge_q, dcfs)):

            if self.rads_list[ndx] != 'H-3':
                # print('max_dilutfac_for_distance_secperm3:', max_dilutfac_for_distance_secperm3)
                inhalation_dose = float(max_dilutfac_for_distance_secperm3) * float(rad) * float(dcf) * float(
                    breathing_rate) * float(1000)

                inhal_dose_per_rad.append(inhalation_dose)
                tot_inhalation_dose += inhalation_dose
            else:
                logging.getLogger("NOTE").info(
                    "Tritium (H-3) in user-defined radionuclide list:: \n{input_desc}".format(
                        input_desc=self.rads_list))
                total_inhalation_dose_tritium = 0
                # for rad, dcf_t in zip(self.discharge_tritium, dcfs):
                inhalation_dose = float(max_dilutfac_for_distance_secperm3) * float(rad) * float(dcf) * float(
                    breathing_rate) * float(1000)
                total_inhalation_dose_tritium += inhalation_dose
                inhal_dose_per_rad.append(total_inhalation_dose_tritium)
                tot_inhalation_dose += total_inhalation_dose_tritium

        return inhal_dose_per_rad

    def ground_shine_dose(self, X1, age=1, max_dilutfac_for_distance_secperm3=None):
        """
            Calculate the ground shine dose for a given release scenario.

            Args:
                X1 (float): Distance from the source in meters.
                age (int, optional): Age of the recipient in years. Default is 1.
                max_dilutfac_for_distance_secperm3 (float, optional): Maximum dilution factor for the given distance in
                    seconds per cubic meter. If not provided, it will be calculated based on the release scenario.

            Returns:
                numpy.ndarray: An array containing the ground shine dose for each radionuclide in the scenario.
        """
        # m3/s
        weathering_corr = self.config['weathering_corr']
        if self.max_dilution_factor is None and self.config['have_met_data']:
            self.file_preprocessing()
            self.met_data_to_tjfd()
            self.missing_correction()
            self.calm_correction_factor_calc()
            self.dilution_per_sector(X1)

        if self.config['have_dilution_factor']:
            max_dilutfac_for_distance_secperm3 = self.dict_max_dilution_factor[X1]
        if not self.config['have_dilution_factor']:
            if max_dilutfac_for_distance_secperm3 is None:
                max_dilutfac_for_distance_secperm3 = self.max_dilution_factor
        # remember: unit is in  Sv-m2/Bq-s
        dcfs = self.dcf_list_ecerman_ground_shine_include_progeny(master_file="library/Dose_ecerman_final.xlsx",
                                                                  sheet_name="surface_dose",
                                                                  age=age,
                                                                  consider_progeny=self.config['consider_progeny'])

        # print('dcfsssssssssssssssssss:', dcfs)
        # deposition velocities of all radionuclides ; see IAEA SRS 19 (page 27)
        list_deposition_vel = self.deposition_velocity_of_rad()

        # lambda
        xls = pd.ExcelFile("library/RadioToxicityMaster.xls")
        name = pd.read_excel(xls, "Inhalation CED Sv per Bq Public")
        name.dropna(axis=0, how='all', inplace=True)
        if self.lambda_of_rads is None:
            self.find_half_life_and_decay_const_radionuclides()

        # applicable only if radionuclide is ['I','Cl',Tc','Cs','Sr']
        # unit is in second dimension
        effective_lambda_frac = self.apply_weathering_correction_gs()

        tot_gs_dose = 0
        gs_dose_per_rad = []

        # for single_plume discharge is in Bq unit (Time integrated)
        if self.config['single_plume']:
            discharge_q = self.config['instantaneous_release_bq_list']
        # for continuous release, discharge in Bq/year unit
        else:
            discharge_q = self.annual_discharge_bq_rad_list

        # dcfs has following form [(corrected_rad_1, uncorrected_rad_1),(corrected_rad_2, uncorrected_rad_2)]
        if self.config['consider_progeny']:
            # take corrected DCF
            for bdx, (discharge_rad_per_year, dcf, effective_lambda) in enumerate(zip(discharge_q, dcfs,
                                                                                      effective_lambda_frac)):
                # DCF is in Sv-m^2/Bq-second
                # deposition rate unit: Bq/m2-second
                # deposition velocity unit: m/s
                # s/m3 * bq/s * m/s
                deposition_rate = float(max_dilutfac_for_distance_secperm3) * (
                        float(discharge_rad_per_year) / (365 * 24 * 3600)) \
                                  * float(list_deposition_vel[bdx])
                # unit Bq/m2
                conc_rad_ground = deposition_rate * float(effective_lambda)
                gs_dose = conc_rad_ground * dcf[0] * float(1000)
                # convert dose to year
                gs_dose = gs_dose * (365 * 24 * 3600)
                logging.getLogger("Ground shine Dose").info(
                    'DCF used for radionuclide {} is {} Sv-m^2/Bq-second [progeny corrected (if required), uncorrected]:'.format(self.rads_list[bdx],
                                                                                          dcf))
                print("Ground shine Dose: DCF used for radionuclide {} is {} Sv-m^2/Bq-second:".format(self.rads_list[bdx],dcf))
                print('Deposition rate for radionuclide {} is {} Bq/m2-second'.format(self.rads_list[bdx],
                                                                                      deposition_rate))
                print('Concentration of radionuclide {} on the ground is {} Bq/m^2:'.format(self.rads_list[bdx],
                                                                                            conc_rad_ground))
                gs_dose_per_rad.append(gs_dose)
                tot_gs_dose += gs_dose

        else:
            # take uncorrected DCF
            for bdx, (discharge_rad_per_year, dcf, effective_lambda) in enumerate(zip(discharge_q, dcfs,
                                                                                      effective_lambda_frac)):
                #    # long-term release
                # DCF is in Sv-m^2/Bq-second
                # deposition rate unit: Bq/m2-second
                # deposition velocity unit: m/s
                # s/m3 * bq/s * m/s
                deposition_rate = float(max_dilutfac_for_distance_secperm3) * (
                        float(discharge_rad_per_year) / (365 * 24 * 3600)) \
                                  * float(list_deposition_vel[bdx])
                # unit Bq/m2
                conc_rad_ground = deposition_rate * float(effective_lambda)
                # multiply by 1000 to convert to mSv from Sv
                gs_dose = conc_rad_ground * dcf[1] * float(1000)
                # convert dose to year
                gs_dose = gs_dose * (365 * 24 * 3600)
                logging.getLogger("Ground shine Dose").info(
                    'DCF used for radionuclide {} is {} Sv-m^2/Bq-second [progeny uncorrected values]:'.format(self.rads_list[bdx],
                                                                                          dcf[0]))
                logging.getLogger("Ground shine Dose").info(
                    'Deposition rate for radionuclide {} is {} Bq/m2-second'.format(self.rads_list[bdx],
                                                                                      deposition_rate))
                logging.getLogger("Ground shine Dose").info(
                    'Concentration of radionuclide {} on the ground is {} Bq/m^2:'.format(self.rads_list[bdx],
                                                                                            conc_rad_ground))

                print("Ground shine Dose: DCF used for radionuclide {} is {} Sv-m^2/Bq-second:".format(
                    self.rads_list[bdx], dcf))
                print('Deposition rate for radionuclide {} is {} Bq/m^2-second'.format(self.rads_list[bdx],
                                                                                      deposition_rate))
                print('Concentration of radionuclide {} on the ground is {} Bq/m^2:'.format(self.rads_list[bdx],
                                                                                            conc_rad_ground))
                # long-term release
                gs_dose_per_rad.append(gs_dose)
                tot_gs_dose += gs_dose

        gs_dose_per_rad = np.array(gs_dose_per_rad)
        self.gs_dose_per_rad = gs_dose_per_rad

        return self.gs_dose_per_rad

    def submersion_dose(self, X1, age=18, max_dilutfac_for_distance_secperm3=None):
        """
            Calculate the submersion dose for a given release scenario.

            Args:
                X1 (float): Distance from the source in meters.
                age (int, optional): Age of the recipient in years. Default is 18.
                max_dilutfac_for_distance_secperm3 (float, optional): Maximum dilution factor for the given distance in
                    seconds per cubic meter. If not provided, it will be calculated based on the release scenario.

            Returns:
                numpy.ndarray: An array containing the submersion dose for each radionuclide in the scenario.
        """

        if self.max_dilution_factor is None and self.config['have_met_data']:
            self.file_preprocessing()
            self.met_data_to_tjfd()
            self.missing_correction()
            self.calm_correction_factor_calc()
            self.dilution_per_sector(X1)

        if self.config['have_dilution_factor']:
            max_dilutfac_for_distance_secperm3 = self.dict_max_dilution_factor[X1]
        if not self.config['have_dilution_factor']:
            if max_dilutfac_for_distance_secperm3 is None:
                max_dilutfac_for_distance_secperm3 = self.max_dilution_factor
        # unit is Sv-m3/Bq-s
        dcfs = self.dcf_list_ecerman_submersion_include_progeny(master_file='library/Dose_ecerman_final.xlsx',
                                                                sheet_name='submersion_dose', age=age)

        # for single_plume discharge is in Bq unit (Time integrated)
        if self.config['single_plume']:
            discharge_q = self.config['instantaneous_release_bq_list']
        # for continuous release, discharge in Bq/year unit
        else:
            discharge_q = self.annual_discharge_bq_rad_list

        tot_submersion_dose = 0
        sub_dose_per_rad = []
        # dcfs has following form [(corrected_rad_1, uncorrected_rad_1),(corrected_rad_2, uncorrected_rad_2)]
        if self.config['consider_progeny']:
            # take corrected DCF
            for bdx, (discharge_rad, dcf) in enumerate(zip(discharge_q, dcfs)):
                logging.getLogger("Submersion Dose").info(
                    'DCF used for radionuclide {} is {} Sv-m3/Bq-s:'.format(self.rads_list[bdx],
                                                                                          dcf))
                print("Submersion Dose: DCF used for radionuclide {} is {} Sv-m3/Bq-second [progeny corrected (if required), uncorrected]:".format(
                    self.rads_list[bdx], dcf))

                if self.config['single_plume']:
                    # converting DCF /sec to /year; for single plume do not convert dcf unit to year
                    # bq-s/m3 * sv/s * m3/bq
                    dcf = dcf[0]
                    submersion_dose = float(max_dilutfac_for_distance_secperm3) * float(discharge_rad) * float(dcf) * \
                                      float(1000)
                else:
                    dcf = dcf[0]

                    # multiplied by 1000 to convert into mSv; s/m3 * Bq/y * Sv-m3/Bq-s = Sv/y * 1000 = mSv/y
                    submersion_dose = float(max_dilutfac_for_distance_secperm3) * float(discharge_rad) * float(dcf) * \
                                      float(1000)
                tot_submersion_dose += submersion_dose
                sub_dose_per_rad.append(submersion_dose)
        else:
            # take uncorrected DCF
            for bdx, (discharge_rad, dcf) in enumerate(zip(discharge_q, dcfs)):
                logging.getLogger("Submersion Dose").info(
                    'DCF used for radionuclide {} is {} Sv-m3/Bq-s [progeny not considered]:'.format(self.rads_list[bdx],
                                                                                          dcf[0]))
                print("Submersion Dose: DCF used for radionuclide {} is {} Sv-m3/Bq-second:".format(
                    self.rads_list[bdx], dcf))
                if self.config['single_plume']:
                    # converting DCF /sec to /year; for single plume do not convert dcf unit to year
                    # bq-s/m3 * sv/s * m3/bq
                    # / (365 * 24 * 3600) check???
                    dcf = dcf[1]
                    submersion_dose = float(max_dilutfac_for_distance_secperm3) * float(discharge_rad) * float(dcf) * \
                                      float(1000)
                else:

                    # DCF Sv-m3/Bq-sec
                    dcf = dcf[1]
                    # multiplied by 1000 to convert into mSv; s/m3 * Bq/y * Sv-m3/Bq-s = Sv/y * 1000 = mSv/y
                    submersion_dose = float(max_dilutfac_for_distance_secperm3) * float(discharge_rad) * float(dcf) * \
                                      float(1000)

                tot_submersion_dose += submersion_dose
                sub_dose_per_rad.append(submersion_dose)
        self.sub_dose_per_rad = np.array(sub_dose_per_rad)
        return self.sub_dose_per_rad

    def ingestion_dose(self, X1, receiver='adult', max_dilutfac_for_distance_secperm3=None):
        """
            Calculate the ingestion dose for a given release scenario.

            Args:
                X1 (float): Distance from the source in meters.
                receiver (str, optional): Receiver type, either 'adult' or 'infant'. Default is 'adult'.
                max_dilutfac_for_distance_secperm3 (float, optional): Maximum dilution factor for the given distance in
                    seconds per cubic meter. If not provided, it will be calculated based on the release scenario.

            Returns:
                numpy.ndarray: An array containing the ingestion dose for each radionuclide in the scenario.
            Raises:
                ValueError: If the receiver type is not 'adult' or 'infant'.
        """

        if self.max_dilution_factor is None and self.config['have_met_data']:
            self.file_preprocessing()
            self.met_data_to_tjfd()
            self.missing_correction()
            self.calm_correction_factor_calc()
            self.dilution_per_sector(X1)

        if self.config['have_dilution_factor']:
            max_dilutfac_for_distance_secperm3 = self.dict_max_dilution_factor[X1]
        if not self.config['have_dilution_factor']:
            if max_dilutfac_for_distance_secperm3 is None:
                max_dilutfac_for_distance_secperm3 = self.max_dilution_factor
            # else:
            #    dilution_factor_sectorwise = self.dilution_per_sector(X1)
            # print('dilution_factor_sectorwise:', dilution_factor_sectorwise)
            # max_dilutfac_for_distance_secperm3 = max(dilution_factor_sectorwise)
            #    max_dilutfac_for_distance_secperm3 = np.array(dilution_factor_sectorwise).T[0].max()
            # max_dilutfac_for_distance_secperm3 = self.max_dilution_factor

        # CHECK AGE OF RECEIVER
        if receiver == 'adult':
            age = 18
        elif receiver == 'infant':
            age = 1
        else:
            raise ValueError('For ingestion dose calculation, the receiver can be either adult or infant.')

        # ASSIGN INGESTION PARAMETER
        if self.config['inges_param_dict']:
            inges_param_dict = dict(self.config['inges_param_dict'])
        else:
            inges_param_dict = {'alpha_wet_crops': 0.3, 'alpha_dry_forage': 3, 't_e_food_crops': 60,
                                't_e_forage_grass': 30, 't_b': 11000, 't_h_wet_crops': 1, 't_h_animal_pasture': 0,
                                't_h_animal_stored_feed': 90,
                                'C_wi': 0, 'f_p': 0.7, 'alpha': 3, 't_e': 30, 't_m': 1, 't_f': 20, 'q_m': 16,
                                'q_w': 0.06, 'q_f': 1.2, 'q_w_meat': 0.004}

        logging.getLogger("Ingestion Dose").info(
            "Parameters used for Ingestion dose calculations: {inges_param_dict}".format(
                inges_param_dict=inges_param_dict))

        # ASSIGN INGESTION PARAMETER TO ADULT
        if self.config['inges_param_dict_adult']:
            inges_param_dict_adult = dict(self.config['inges_param_dict_adult'])
        else:
            inges_param_dict_adult = {'DID_veg': 76.7, 'DID_milk': 182.5, 'DID_meat': 14.6, 'DID_fish': 18.3,
                                      'DID_water_and_beverage': 0.73}

        # ASSIGN INGESTION PARAMETER TO INFANT
        if self.config['inges_param_dict_infant']:
            inges_param_dict_infant = dict(self.config['inges_param_dict_infant'])
        else:
            inges_param_dict_infant = {'DID_veg': 78.5, 'DID_milk': 146, 'DID_meat': 1.2, 'DID_fish': 1.5,
                                       'DID_water_and_beverage': 0.26}

        # PRINT USED PARAMETERS IN LOG FILE
        logging.getLogger("Ingestion Dose").info(
            "Additional Parameters used for Ingestion dose calculations (Adult): {inges_param_dict_adult}".format(
                inges_param_dict_adult=inges_param_dict_adult))

        logging.getLogger("Ingestion Dose").info(
            "Additional Parameters used for Ingestion dose calculations (Infant): {inges_param_dict_infant}".format(
                inges_param_dict_infant=inges_param_dict_infant))

        # CHECK THE SOURCE TERM
        # for single_plume discharge is in Bq unit (Time integrated)
        if self.config['single_plume']:
            discharge_q = self.config['instantaneous_release_bq_list']
            day_discharge_bq_rad_list = discharge_q
        # for continuous release, discharge in Bq/year unit
        else:
            discharge_q = self.annual_discharge_bq_rad_list
            day_discharge_bq_rad_list = [float(each) / 365 for each in discharge_q]

        # deposition velocity (unit m/s); see IAEA SRS 19 (page 27)
        list_deposition_vel = self.deposition_velocity_of_rad()
        logging.getLogger("Ingestion Dose").info(
            "Used deposition velocities in m/s: {list_deposition_vel}".format(
                list_deposition_vel=list_deposition_vel))
        # (ref. Page 10, SRS 19) The deposition rate d^dot_i is used to calculate the radionuclide concentration on
        # vegetation owing to direct contamination; unit # Bq/m2/d
        deposition_rate_per_day = [d * list_deposition_vel[idx] * max_dilutfac_for_distance_secperm3 for idx, d in
                                   enumerate(day_discharge_bq_rad_list)]
        print('deposition_rate_per_day for the given list of radionclides:', deposition_rate_per_day)

        # GET DECAY CONSTANT; unit /second
        if self.lambda_of_rads is None:
            self.find_half_life_and_decay_const_radionuclides()
        lambda_i = self.lambda_of_rads
        # converted to /day as all variables in calculations in ingestion dose is in /d unit.
        lambda_i = [_ * 24 * 3600 for _ in lambda_i]
        # PERFORM WEATHERING CORRECTION:
        # The radionuclide concentration in vegetation may be reduced by a variety of
        # processes. These include radioactive decay, wash-off of previously intercepted material by rain or
        # irrigation, surface abrasion and leaf bending from the action of the wind, resuspension, tissue ageing,
        # leaf fall or herbivore grazing, addition of new tissue (growth dilution), volatilization or evaporation.
        # Losses other than radioactive decay are normally described using an aggregated parameter in the form of
        # a first order rate constant lw. Table 7 page 63-64, SRS 19 A default value of lw is given in Table VII for
        # estimating the removal of radionuclides from vegetation based on a half-life of 14 days.
        # noble_gas_tritium_c14 = ['H', 'C', 'He', 'Ne', 'Ar', 'Kr', 'Xe', 'Rn']
        # noble_gas = ['He', 'Ne', 'Ar', 'Kr', 'Xe', 'Rn']

        # self.effective_lambda_of_rads = self.ingestion_weathering_correction()

        ############################### VEG ROUTE #####################################################################

        # alpha (m2/kg) = fraction of deposited activity intercepted by the edible portion
        # of vegetation per unit mass ; source: table VIII, page 64 of SRS 19
        alpha_wet_crops = inges_param_dict['alpha_wet_crops']
        alpha_dry_forage = inges_param_dict['alpha_dry_forage']

        # F_v (Bq/kg dry soil) = concentration for uptake of radionuclide from soil by edible parts of crops.
        # source: page 67, Table XI of SRS 19; Take value from Table B.4
        lambda_s_per_d_list, lambda_w_per_d_list, f_v1_list, f_v2_list, Fm_Milk_d_per_L_list, Ff_Meat_d_per_kg_list, notransfer_factor_rad = \
            self.fv_list_ecerman_ingestion(master_file="library/Dose_ecerman_final.xlsx", sheet_name="eco_param")
        # C_vi1 (Bq/Kg) = vegetation consumed by grazing animal or human; other than root uptake

        # THIS IS WEATHERING CORRECTION FOR INGESTION USING ECO-PARAM ####
        # converted to /day from /second as lambda_w_per_d_list and lambda_s_per_d_list has unit /day
        # lambda_eiv_list and lambda_eis_list have unit /day
        lambda_eiv_list, lambda_eis_list = self.ingestion_weathering_correction_real(lambda_w_per_d_list,
                                                                                     lambda_s_per_d_list)
        print('lambda_w_per_d_list,lambda_s_per_d_list:', lambda_w_per_d_list, lambda_s_per_d_list,
              lambda_eiv_list, lambda_eis_list)

        # rho = standardised surface density for the effective root zone in soil (kg/m2 dry weight soil)
        # from TABLE IX. EFFECTIVE SURFACE SOIL DENSITY FOR SCREENING PURPOSES (IAEA SRS 19)
        rho_pasture_depth_lt_11, rho_crop_depth_ge_11 = self.effective_surface_soil_density_rho()

        # unit: days #source: table VIII, page 64 of SRS 19
        t_e_food_crops = inges_param_dict['t_e_food_crops']

        t_b = inges_param_dict['t_b']  # t_b = duration of discharge of material; for 30 years it is 11000 days.
        # t_h_wet_crops: unit => day; delay time (time interval between harvest and consumption of food; default 1 day)
        t_h_wet_crops = inges_param_dict['t_h_wet_crops']

        # FINDING TOTAL CONCENTRATION OF RADIONUCLIDES ON THE VEGETATION AT THE TIME OF CONSUMPTION
        # STEP 1: compute Cvi and store it in Cvi_list. Cvi (Bq/Kg) = Bq/kg dry for vegetation consumed by grazing
        # STEP 2: MULTIPLY Cvi with np.exp(-lambda_i[ndx] * t_h_wet_crops (default value = 1 day)
        # animals; bq/kg fresh matter for vegetation consumed by humans.
        Cvi_list = []
        # deposition_rate_per_day = [1]
        for ndx, each in enumerate(self.element_list):
            if each not in ['H', 'C']:
                # STEP 1
                # C_vi1: due to direct contamination of nuclide in and on vegetation; Bq/Kg
                C_vi1 = (deposition_rate_per_day[ndx] * alpha_wet_crops * (
                        1 - np.exp(-lambda_eiv_list[ndx] * t_e_food_crops))) / (lambda_eiv_list[ndx])
                logging.getLogger("Ingestion Dose").info(
                    'Concentration of {} on vegetation surfaces owing to direct '
                    'deposition from the atmosphere is {} Bq/Kg fresh weight in food crops.'.format(each, C_vi1))
                print('Concentration of {} on vegetation surfaces owing to direct '
                      'deposition from the atmosphere is {} Bq/Kg fresh weight in food crops.'.format(each, C_vi1))
                # C_si (Bq/Kg) =  concentration of radionuclide in food crops from uptake from dry soil;
                # rho_crop_depth_ge_11 is taken as default; may be added an option later
                C_si = (deposition_rate_per_day[ndx] * (1 - np.exp(-lambda_eis_list[ndx] * t_b))) / (
                        lambda_eis_list[ndx] * rho_crop_depth_ge_11)
                logging.getLogger("Ingestion Dose").info(
                    'Concentration of radionuclide {} in food crops from uptake from soil '
                    'is {} Bq/Kg dry weight of soil.'.format(each, C_si))
                print('Concentration of radionuclide {} in food crops from uptake from soil '
                      'is {} Bq/Kg dry weight of soil.'.format(each, C_si))
                # C_vi2: Bq/kg dry for vegetation consumed by grazing animals; bq/kg fresh matter for vegetation
                # consumed by humans.
                C_vi2 = f_v2_list[ndx] * C_si
                logging.getLogger("Ingestion Dose").info(
                    'The concentration of radionuclide {} in vegetation resulting from uptake from soil  '
                    'is {} Bq/Kg fresh weight.'.format(each, C_vi2))
                print('The concentration of radionuclide {} in vegetation resulting from uptake from soil  '
                      'is {} Bq/Kg fresh weight.'.format(each, C_vi2))
                # STEP 2
                Cvi = (C_vi1 + C_vi2) * np.exp(-lambda_i[ndx] * t_h_wet_crops)
                logging.getLogger("Ingestion Dose").info(
                    'C_food_crops: The concentration of {} in food crops from both direct deposition and uptake '
                    'from soil is {} Bq/Kg fresh weight.'.format(each, Cvi))
                print('C_food_crops: The concentration of {} in food crops from both direct deposition and uptake '
                      'from soil is {} Bq/Kg fresh weight.'.format(each, Cvi))
                Cvi_list.append(Cvi)

        ###################  MILK AND MEAT ROUTE ##################################################################
        # conc in animal feed
        t_h_animal_pasture = inges_param_dict['t_h_animal_pasture']
        t_h_animal_stored_feed = inges_param_dict['t_h_animal_stored_feed']
        # f_p: fraction of the year that animal consumes fresh pasture vegetation
        f_p = inges_param_dict['f_p']
        # c_pi = concentration of radionuclide in stored feeds (Bq/Kg); t_h_delay_time in day
        # t_h_delay_time = 90
        t_e = inges_param_dict['t_e']  # check for correctness
        t_m = inges_param_dict['t_m']  # unit: day# average time between collection and human consumption of milk
        t_f = inges_param_dict[
            't_f']  # average time between collection and human consumption of meat ; 20 d (IAEA SRS 19)
        q_m = inges_param_dict['q_m']  # amount of feed (in dry matter) consumed per day (kg/d)
        q_w = inges_param_dict['q_w']  # amount of water consumed by animal per day (m3/d) #Table B1 B2 page 66 (ECPDA)
        C_wi = inges_param_dict['C_wi']  # concentration of radionuclide in water (in Bq/m3) #check
        q_f = inges_param_dict[
            'q_f']  # amount of feed consumed by animal (Kg/d) ; goat and sheep ; meat producing animal
        q_w_meat = inges_param_dict['q_w_meat']  # water intake of meat producing animal (m3/d)
        # t_e = the period that the crops are exposed to contamination during growing season. unit: day
        t_e_forage_grass = inges_param_dict['t_e_forage_grass']
        # concentration of radionuclide in stored feeds (Bq/Kg, dry matter)
        Cpi_list = []
        # concentration of radionuclide in animal feed (Bq/Kg, dry matter)
        Cai_list = []
        # concentration of radionuclide in milk (Bq/l)
        Cmi_list = []
        # concentration of radionuclide in flesh (Bq/Kg)
        Cfi_list = []

        for ndx, each in enumerate(self.element_list):
            if each not in ['H', 'C']:
                # Bq/kg DW: Cvi calculation for dry_forage (not wet!!!) for entering into meat and milk route
                C_vi1 = (deposition_rate_per_day[ndx] * alpha_dry_forage * (
                        1 - np.exp(-lambda_eiv_list[ndx] * t_e_forage_grass))) / (lambda_eiv_list[ndx])
                logging.getLogger("Ingestion Dose").info(
                    'The concentration of {} on forage owing to direct deposition is {} Bq/Kg dry weight forage.'
                    ''.format(each, C_vi1))
                print('The concentration of {} on forage owing to direct deposition is {} Bq/Kg dry weight forage.'
                      ''.format(each, C_vi1))
                # Bq/kg
                C_si = (deposition_rate_per_day[ndx] * (1 - np.exp(-lambda_eis_list[ndx] * t_b))) / (
                        lambda_eis_list[ndx] * rho_pasture_depth_lt_11)
                logging.getLogger("Ingestion Dose").info(
                    'The concentration of {} is pasture soil is {} Bq/Kg dry weight of soil.'
                    ''.format(each, C_si))
                print('The concentration of {} is pasture soil is {} Bq/Kg dry weight of soil.'
                      ''.format(each, C_si))
                # Bq/Kg
                C_vi2 = f_v1_list[ndx] * C_si
                logging.getLogger("Ingestion Dose").info(
                    'The concentration of {} through uptake from soil by pasture vegetation is {}'.format(each, C_vi2))
                print(
                    'The concentration of {} through uptake from soil by pasture vegetation is {}'.format(each, C_vi2))
                # Bq/Kg
                Cvi_animal = (C_vi1 + C_vi2) * np.exp(-lambda_i[ndx] * t_h_animal_pasture)
                logging.getLogger("Ingestion Dose").info(
                    'C_pasture: The total contamination of {} in pasture from direct deposition and uptake from soil '
                    'is {} Bq/Kg dry weight forage. '
                    ''.format(each, Cvi_animal))
                print(
                    'C_pasture: The total contamination of {} in pasture from direct deposition and uptake from soil '
                    'is {} Bq/Kg dry weight forage. '
                    ''.format(each, Cvi_animal))

                # Concentrations in stored feed and average concentrations for feeds
                Cpi = (C_vi1 + C_vi2) * np.exp(-lambda_i[ndx] * t_h_animal_stored_feed)
                logging.getLogger("Ingestion Dose").info(
                    'The concentration of {} with storage time of {} day '
                    'on stored feed is {} Bq/Kg dry weight forage'.format(each, t_h_animal_stored_feed, Cpi))
                print('The concentration of {} with storage time of {} day '
                      'on stored feed is {} Bq/Kg dry weight forage'.format(each, t_h_animal_stored_feed, Cpi))
                Cpi_list.append(Cpi)
                # concentration of radionuclide in animal feed (Bq/Kg, dry matter)
                forage_feed_percent = f_p * 100
                C_ai = f_p * Cvi_animal + (1 - f_p) * Cpi
                logging.getLogger("Ingestion Dose").info(
                    'Assuming that both dairy and beef cattle are fed a diet of {} fresh forage '
                    'the average annual concentration of {} in animal feed is {} Bq/kg dry weight.'
                    ''.format(each, forage_feed_percent, C_ai))
                print('Assuming that both dairy and beef cattle are fed a diet of {} fresh forage '
                      'the average annual concentration of {} in animal feed is {} Bq/kg dry weight.'
                      ''.format(each, forage_feed_percent, C_ai))
                Cai_list.append(C_ai)

                # milk route; # concentration of radionuclide in milk (Bq/l)
                C_mi = Fm_Milk_d_per_L_list[ndx] * ((C_ai * q_m) + (C_wi * q_w)) * np.exp(-lambda_i[ndx] * t_m)
                logging.getLogger("Ingestion Dose").info(
                    'C_milk: The concentration of {} in milk is {} Bq/L.'.format(each, C_mi))
                print('C_milk: The concentration of {} in milk is {} Bq/L.'.format(each, C_mi))
                if C_wi == 0 or q_w == 0:
                    print('Note: It is assumed that the concentration of radionuclide in water is negligible.')
                Cmi_list.append(C_mi)

                # meat route
                C_fi = Ff_Meat_d_per_kg_list[ndx] * ((C_ai * q_f) + (C_wi * q_w_meat)) * np.exp(-lambda_i[ndx] * t_f)
                logging.getLogger("Ingestion Dose").info(
                    'C_meat: The concentration of {} in meat is {} Bq/kg.'.format(each, C_fi))
                print('C_meat: The concentration of {} in meat is {} Bq/kg.'.format(each, C_fi))
                if q_f == 12:
                    print("Note: Large animals are considered as meat source.")
                if q_f == 1.2:
                    print("Note: Small animals are considered as meat source.")
                Cfi_list.append(C_fi)

                # fish route (TODO)
        # loading excel for dcf extraction due to ingestion
        dcfs = self.dcf_list_ingestion(master_file="library/Dose_ecerman_final.xlsx", sheet_name="ingestion_gsr3",
                                       age=age)

        if self.config['long_term_release']:
            logging.getLogger("Ingestion Dose").info(
                "For long-term release:: \n{input_desc}".format(input_desc=self.config['long_term_release']))
            # dietary intake data
            if receiver == 'adult':
                # convert to year only for long-term release;
                DID_veg = inges_param_dict_adult['DID_veg'] * 365  # 383.3 #default but ask to user; kg/year
                DID_milk = inges_param_dict_adult['DID_milk'] * 365  # l/year
                DID_meat = inges_param_dict_adult['DID_meat'] * 365  # kg/year
                DID_fish = inges_param_dict_adult['DID_fish'] * 365  # kg/year
                DID_water_and_beverage = inges_param_dict_adult['DID_water_and_beverage'] * 365  # m3/year
                # TO DO: can be directly read from Dose_ecerman_final file which provides values for other age group
                # as well dcf_hto = 1.8E-11  # unit Sv/Bq taken from "library/Dose_ecerman_final.xlsx",
                # sheet_name ="ingestion_gsr3" dcf_obt = 4.2E-11  # unit Sv/Bq taken from
                # "library/Dose_ecerman_final.xlsx", sheet_name= "ingestion_gsr3"
            elif receiver == 'infant':
                DID_veg = inges_param_dict_infant['DID_veg'] * 365  # default but ask to user; kg/year
                DID_milk = inges_param_dict_infant['DID_milk'] * 365  # l/year
                DID_meat = inges_param_dict_infant['DID_meat'] * 365  # kg/year
                DID_fish = inges_param_dict_infant['DID_fish'] * 365  # kg/year
                DID_water_and_beverage = inges_param_dict_infant['DID_water_and_beverage'] * 365  # m3/year
                # dcf_hto = 6.4E-11  # unit Sv/Bq taken from "library/Dose_ecerman_final.xlsx", sheet_name="ingestion_gsr3"
                # dcf_obt = 1.2E-10  # unit Sv/Bq taken from "library/Dose_ecerman_final.xlsx", sheet_name="ingestion_gsr3"
            else:
                raise ValueError('For ingestion dose calculation, the receiver can be either adult or infant.')

        if self.config['single_plume']:
            consumption_time_frac = 1
            # dietary intake data
            if receiver == 'adult':
                DID_veg = inges_param_dict_adult[
                              'DID_veg'] * consumption_time_frac  # 383.3 #default but ask to user; kg/year
                DID_milk = inges_param_dict_adult['DID_milk'] * consumption_time_frac  # l/year
                DID_meat = inges_param_dict_adult['DID_meat'] * consumption_time_frac  # kg/year
                DID_fish = inges_param_dict_adult['DID_fish'] * consumption_time_frac  # kg/year
                DID_water_and_beverage = inges_param_dict_adult[
                                             'DID_water_and_beverage'] * consumption_time_frac  # m3/year

            elif receiver == 'infant':
                DID_veg = inges_param_dict_infant['DID_veg'] * consumption_time_frac  # default but ask to user; kg/year
                DID_milk = inges_param_dict_infant['DID_milk'] * consumption_time_frac  # l/year
                DID_meat = inges_param_dict_infant['DID_meat'] * consumption_time_frac  # kg/year
                DID_fish = inges_param_dict_infant['DID_fish'] * consumption_time_frac  # kg/year
                DID_water_and_beverage = inges_param_dict_infant[
                                             'DID_water_and_beverage'] * consumption_time_frac  # m3/year

            else:
                raise ValueError('For ingestion dose calculation, the receiver can be either adult or infant.')

        # THE CURIOUS CASE OF TRITIUM ##################################################################################
        # concentration of tritium in terrestrial plant products
        def conc_tritium_in_terrestrial_plant(max_dilution_factor, C_HTO_atm, climate='Continental',
                                              veg_type='leafy_vegetables',
                                              CR_s=0.23, gamma=0.909, R_p=0.54):
            # C_HTO_atm = Bq/year
            # climate : Mediterranean, Continental, Maritime, Arctic (for relative and absolute humidity values)
            # C_tfwt is the tritium concentration in the tissue free leaf water (Bq/m3)
            # gamma = ratio of vapour pressure of tritiated water (HTO) to normal water (H2O). The value is 0.909.
            # C_am = annual average tritium concentration in air moisture.
            # C_HTO_atm is the annual average concentration of tritium in ground level air estimated from release rate
            # and atmospheric dispersion data in (Bq/m3)
            # H_a is the average absolute humidity (Kg/m3)
            # C_HTO_atm = WCp * C_tfwt; modeling or observed data

            # e_s = saturation vapour pressure in air (Pa)
            # T_a = air temperature
            # RH = relative humidity
            # H_a  = (0.00217 * e_s * RH)/T_a

            # converted to Bq/second
            C_HTO_atm = (C_HTO_atm * max_dilution_factor) / (365 * 24 * 3600)

            climate_humidity = {'Mediterranean': [34, 0.0115, 0.6], 'Continental': [48, 0.0087, 0.71],
                                'Maritime': [50, 0.0078, 0.795], 'Arctic': [50, 0.0067, 0.73]}

            # for all others, WCp is the average value from table 3, IAEA tecdoc 1616.
            veg_type_data = {'leafy_vegetables': [0.51, 0.92], 'non_leafy_vegetables': [0.53, 0.92],
                             'root_crops': [0.52, 0.87], 'all_others': [0.56, 0.495]}

            latitude, H_a, RH = climate_humidity[climate]

            WEQ, WCp = veg_type_data[veg_type]

            C_am = C_HTO_atm / H_a

            # CR_s is the ratio of concentrations of tritium in soil water to that in air moisture (unitless)
            # CR_s =  2.3 * 10-1; Geometric mean taken from IAEA Tecdoc 1616
            C_sw = CR_s * C_am
            # Bq/L
            C_tfwt = ((RH * C_am) + (1 - RH) * C_sw) / gamma

            # C_pfw_obt = the concentration of OBT in Bq/Kg fw of plant product.
            C_pfw_obt = (1 - WCp) * WEQ * R_p * C_tfwt

            # WEQ is the water equivalent factor, Litre of water produced per Kg of dry matter combusted (L/Kg)
            # WCp = fractional water content in the biota
            # R_p = ratio of concetration of OBT (Bq/L) to that of TFWT (Bq/L) in the biota.

            #  C_pfw_hto: The HTO concentration in the fresh weight (FW) plant (Bq kg-1 FW)
            C_pfw_hto = WCp * C_tfwt

            return WEQ, WCp, C_tfwt, C_pfw_hto, C_pfw_obt

        # tritium concentration in terrestrial animal products
        def conc_tritium_in_terrestrial_animal(C_tfwt, C_f_HTO=2.5E5, animal_product='cow_milk',
                                               animal_feed_type='leafy_vegetables',
                                               R_p=0.54):
            """
                Calculate the concentration of tritium in terrestrial plants.

                Args:
                    max_dilution_factor (float): Maximum dilution factor for the given distance.
                    C_HTO_atm (float): Annual average concentration of tritium in ground-level air (Bq/m3).
                    climate (str, optional): Climate type, one of 'Mediterranean', 'Continental', 'Maritime', or 'Arctic'.
                        Defaults to 'Continental'.
                    veg_type (str, optional): Type of vegetation, one of 'leafy_vegetables', 'non_leafy_vegetables',
                        'root_crops', or 'all_others'. Defaults to 'leafy_vegetables'.
                    CR_s (float, optional): Ratio of concentrations of tritium in soil water to that in air moisture.
                        Defaults to 0.23.
                    gamma (float, optional): Ratio of vapour pressure of tritiated water (HTO) to normal water (H2O).
                        Defaults to 0.909.
                    R_p (float, optional): Ratio of concentration of OBT (Bq/L) to that of TFWT (Bq/L) in the biota.
                        Defaults to 0.54.

                Returns:
                    tuple: A tuple containing the following:
                        - WEQ (float): Water equivalent factor, Litre of water produced per Kg of dry matter combusted (L/Kg).
                        - WCp (float): Fractional water content in the biota.
                        - C_tfwt (float): Tritium concentration in the tissue free leaf water (Bq/m3).
                        - C_pfw_hto (float): HTO concentration in the fresh weight (FW) plant (Bq/kg FW).
                        - C_pfw_obt (float): Concentration of OBT in Bq/kg fresh weight of plant product.
            """

            # C_f_HTO = average HTO concentration in ingested water (Bq/L); it can be computed from aquatic modeling

            # C_TFWT can be obtained using conc_tritium_in_terrestrial_plant function
            # animal_product has following options: cow_milk, goat_milk, goat_meat, lamb_meat, beef_meat, pork_meat, broiler_meat, egg
            # C_afw_T_HTO is the total tritium concentration in the animal product from HTO intake (Bq.Kg-1 FW)
            # CR_a_HTO is the concentration ratio for HTO intake (Bq/Kg FW/Bq.L-1)
            # AT  THE MOMENT 2.5E5 is taken (Bq/m3); no validation

            CR_a_HTO_data = {'cow_milk': [0.87], 'goat_milk': [0.83],
                             'goat_meat': [0.67], 'lamb_meat': [0.78], 'beef_meat': [0.66],
                             'pork_meat': [0.67], 'broiler_meat': [0.76], 'egg': [0.66]}

            CR_a_HTO = CR_a_HTO_data[animal_product][0]

            CR_a_OBT_data = {'cow_milk': [0.24], 'goat_milk': [0.32],
                             'goat_meat': [0.43], 'lamb_meat': [0.55], 'beef_meat': [0.40],
                             'pork_meat': [0.64], 'broiler_meat': [0.5], 'egg': [0.64]}

            CR_a_OBT = CR_a_OBT_data[animal_product][0]

            # for all others, WCp is the average value from table 3, IAEA tecdoc 1616.
            animal_feed_type_data = {'leafy_vegetables': [0.51, 0.92], 'non_leafy_vegetables': [0.53, 0.92],
                                     'root_crops': [0.52, 0.87], 'all_others': [0.56, 0.495]}

            WEQ, WCp = animal_feed_type_data[animal_feed_type]

            C_afw_T_HTO = CR_a_HTO * C_f_HTO

            # C_f_OBT = average OBT concentration in feed (Bq Kg-1 DW)
            # WEQ = water equivalent factor, litre of water produced per kg of dry matter combusted (L/Kg)
            # R_p = ratio of concentration of OBT (Bq.L-1) to that of TFWT (Bq.L-1) in the biota
            # C_TFWT is the tritium concentration in the tissue free leaf water (Bq.L-1)
            C_f_OBT = WEQ * R_p * C_tfwt

            # C_afw_T_OBT = The total tritium concentration in the animal product from OBT intake (Bq/Kg FW)
            # CR_a_OBT = concentration ratio for OBT intake (Bq.Kg-1 FW / Bq.Kg-1 DW)
            # C_f_OBT = average OBT concentration in feed (Bq Kg-1 DW); is a weighted
            # average that includes uncontaminated as well as contaminated feed since local sources supply
            # only a fraction of the total animal feed in modern industrial farming
            C_afw_T_OBT = CR_a_OBT * C_f_OBT

            return C_afw_T_HTO, C_f_OBT, C_afw_T_OBT

        # C-14 concentration in terrestrial plants
        def conc_c14_in_terrestrial_plants(C_air, veg_type='leafy_vegetables', S_air=0.20):
            """
                Calculate the concentration of C-14 in terrestrial plants.

                Args:
                    C_air (float): Concentration of C-14 in air (Bq/m3).
                    veg_type (str, optional): Type of vegetation, one of 'leafy_vegetables', 'non_leafy_vegetables',
                        'root_crops', or 'all_others'. Defaults to 'leafy_vegetables'.
                    S_air (float, optional): Concentration of stable carbon in air (gC/m3). Defaults to 0.20.

                Returns:
                    float: Concentration of C-14 in plant (Bq/kg FW).
            """
            # C_pfw is the concentration of C-14 in plant  (Bq/Kg FW)
            # S_p = the concentration of stable carbon in the plant (gC/Kg FW)
            # C_air = concentration of C-14 in air (Bq/m3)
            # S_air  = concentration of stable carbon in air (gC/m3), presently about 0.20 g/m3
            # for all others, S_p is the average value from table , IAEA tecdoc 1616.
            veg_type_S_p_data = {'leafy_vegetables': [30], 'non_leafy_vegetables': [30],
                                 'root_crops': [46], 'all_others': [219]}
            S_p = veg_type_S_p_data[veg_type][0]
            C_pfw = (C_air * S_p) / S_air
            return C_pfw

        # C-14 concentration in terrestrial animal products (Bq/Kg FW)
        # C_pfw = conc_c14_in_terrestrial_plants(C_air, veg_type = 'leafy_vegetables', S_air=0.20)
        def conc_c14_in_terrestrial_animal(C_pfw, animal_feed_type='leafy_vegetables', animal_product='cow_milk',
                                           f_c=1):
            """
                Calculate the concentration of C-14 in terrestrial animal products.

                Args:
                    C_pfw (float): Concentration of C-14 in plant (Bq/kg FW).
                    animal_feed_type (str, optional): Type of animal feed, one of 'leafy_vegetables', 'non_leafy_vegetables',
                        'root_crops', or 'all_others'. Defaults to 'leafy_vegetables'.
                    animal_product (str, optional): Type of animal product, one of 'cow_milk', 'goat_milk', 'goat_meat',
                        'lamb_meat', 'beef_meat', 'pork_meat', 'broiler_meat', or 'egg'. Defaults to 'cow_milk'.
                    f_c (float, optional): Fraction of animal feed that is contaminated. Defaults to 1.

                Returns:
                    float: Concentration of C-14 in animal product (Bq/kg FW).
            """

            # f_c = fraction of animal feed that is contaminated; conservatively taken 1
            # S_a = concentration of stable carbon in the animal product (gC/Kg FW)
            # S_p = the concentration of stable carbon in the plant (gC/Kg FW)
            # C_air = concentration of C-14 in air (Bq/m3)
            # S_air  = concentration of stable carbon in air (gC/m3), presently about 0.20 g/m3
            animal_feed_type_S_p_data = {'leafy_vegetables': [30], 'non_leafy_vegetables': [30],
                                         'root_crops': [46], 'all_others': [219]}
            S_p = animal_feed_type_S_p_data[animal_feed_type][0]

            # Table 12, IAEA Tecdoc 1616
            S_a_data = {'cow_milk': [65], 'goat_milk': [71],
                        'goat_meat': [170], 'lamb_meat': [280], 'beef_meat': [200],
                        'pork_meat': [300], 'broiler_meat': [150], 'egg': [160]}

            S_a = S_a_data[animal_product][0]

            C_apf = (f_c * C_pfw * S_a) / S_p
            return C_apf

        # C-14 concentration in fish (Bq/Kg FW) # Not implemented yet # TODO
        def conc_c14_in_fish(C_air, S_f=120):
            """
                Calculate the concentration of C-14 in fish.

                Args:
                    C_air (float): Concentration of C-14 in air (Bq/m3).
                    S_f (float, optional): Concentration of stable carbon in the fish (gC/kg FW). Defaults to 120.

                Returns:
                    float: Concentration of C-14 in fish (Bq/kg FW).
            """

            # C_DIC = C14 concentration in dissolved inorganic carbon in water column
            # S_f = conc of stable carbon in the fish.
            C_ffw = C_DIC * S_f
            return C_ffw

        ingestion_dose = []
        for ndx, each in enumerate(self.rads_list):
            if each not in ['H-3', 'C-14']:
                # dose due to consumption of leafy vegetables #mSv/year;
                dose_veg = Cvi_list[ndx] * dcfs[ndx] * DID_veg * float(1000)
                # dose due to consumption of milk #mSv/year
                dose_milk = Cmi_list[ndx] * dcfs[ndx] * DID_milk * float(1000)
                # dose due to consumption of meat #mSv/year
                dose_meat = Cfi_list[ndx] * dcfs[ndx] * DID_meat * float(1000)
                # print('ings:', dose_veg, dose_milk, dose_meat)
                ingestion_dose.append((dose_veg, dose_milk, dose_meat))

            # self.ingestion_dose_vals = np.array(ingestion_dose)
            print('dcfs_ing:', dcfs)
            if each == 'H-3':
                logging.getLogger("NOTE: Ingestion").info(
                    "Tritium (H-3) in user-defined radionuclide list:: \n{input_desc}".format(
                        input_desc=self.rads_list))
                print('dcfs:', dcfs, list(dcfs))
                if each == 'H-3' and len(self.rads_list) == 1:
                    dcf_hto = dcfs[0][0]
                    dcf_obt = dcfs[0][1]
                else:
                    dcf_hto = [_ for _ in list(dcfs) if isinstance(_, list)][0][0]
                    dcf_obt = [_ for _ in list(dcfs) if isinstance(_, list)][0][1]
                if self.config['long_term_release']:
                    atm_discharge_HTO = [self.config['annual_discharge_bq_rad_list'][ndx] for ndx, _ in
                                         enumerate(self.rads_list) if _ == 'H-3'][0]
                if self.config['single_plume']:
                    atm_discharge_HTO = [self.config['instantaneous_release_bq_list'][ndx] for ndx, _ in
                                         enumerate(self.rads_list) if _ == 'H-3'][0]

                # for veg
                ingestion_dose_tritium = []
                sum_dose_veg_tritium = 0
                C_HTO_atm = atm_discharge_HTO
                for veg_type in [self.config['veg_type_list']]:
                    WEQ, WCp, C_tfwt, C_pfw_hto, C_pfw_obt = conc_tritium_in_terrestrial_plant(
                        max_dilutfac_for_distance_secperm3, C_HTO_atm, climate=self.config['climate'],
                        veg_type=veg_type, CR_s=0.23, gamma=0.909, R_p=0.54)

                    # to convert it in mSv/y
                    dose_veg_tritium_hto = C_pfw_hto * dcf_hto * DID_veg * float(1000)
                    dose_veg_tritium_obt = C_pfw_obt * dcf_obt * DID_veg * float(1000)
                    sum_hto_obt = dose_veg_tritium_hto + dose_veg_tritium_obt
                    sum_dose_veg_tritium += sum_hto_obt

                # for animal (milk or meat); fish and egg is not yet implemented
                sum_dose_milk_tritium = 0
                sum_dose_meat_tritium = 0
                for animal_prod in self.config['animal_product_list_for_tritium']:
                    # get C_tfwt based on animal feed type (one of the input in input.yaml file)
                    WEQ, WCp, C_tfwt, C_pfw_hto, C_pfw_obt = conc_tritium_in_terrestrial_plant(
                        max_dilutfac_for_distance_secperm3,
                        C_HTO_atm=C_HTO_atm,
                        climate=self.config['climate'],
                        veg_type=self.config['animal_feed_type'],
                        CR_s=0.23, gamma=0.909,
                        R_p=0.54)

                    # use C_tfwt and get concentration of tritium for and in terrestrial animals
                    # note: C_f_HTO is user-defined parameter (from observation or aquatic modeling)
                    C_afw_T_HTO, C_f_OBT, C_afw_T_OBT = conc_tritium_in_terrestrial_animal(C_tfwt=C_tfwt, C_f_HTO=2.5E5,
                                                                                           animal_product=animal_prod,
                                                                                           animal_feed_type=self.config[
                                                                                               'animal_feed_type'],
                                                                                           R_p=0.54)
                    # concentration to dose using dcf for hto and obt
                    if 'DID_milk' in inges_param_dict_adult and any(
                            x in ['cow_milk', 'goat_milk'] for x in self.config['animal_product_list_for_tritium']):
                        # to convert it in mSv/y
                        dose_terres_anim_tritium_hto = C_afw_T_HTO * max_dilutfac_for_distance_secperm3 \
                                                       * dcf_hto * DID_milk * float(1000)
                        dose_terres_anim_tritium_obt = C_afw_T_OBT * max_dilutfac_for_distance_secperm3 \
                                                       * dcf_obt * DID_milk * float(1000)
                        sum_hto_obt_animal_milk_tritium = dose_terres_anim_tritium_hto + dose_terres_anim_tritium_obt
                        sum_dose_milk_tritium += sum_hto_obt_animal_milk_tritium

                    if 'DID_meat' in inges_param_dict_adult and any(
                            x in ['cow_meat', 'goat_meat', 'lamb_meat', 'beef_meat', 'broiler_meat', 'pork_meat'] for x
                            in self.config['animal_product_list_for_tritium']):
                        # to convert it in mSv/y
                        dose_terres_anim_tritium_hto = C_afw_T_HTO * max_dilutfac_for_distance_secperm3 \
                                                       * dcf_hto * DID_meat * float(1000)
                        dose_terres_anim_tritium_obt = C_afw_T_OBT * max_dilutfac_for_distance_secperm3 \
                                                       * dcf_obt * DID_meat * float(1000)
                        sum_hto_obt_animal_meat_tritium = dose_terres_anim_tritium_hto + dose_terres_anim_tritium_obt
                        sum_dose_meat_tritium += sum_hto_obt_animal_meat_tritium

                ingestion_dose_tritium.append(
                    ([sum_dose_veg_tritium], [sum_dose_milk_tritium], [sum_dose_meat_tritium]))

            if each == 'C-14':
                logging.getLogger("NOTE: Ingestion").info(
                    "Carbon-14 (C-14) in user-defined radionuclide list:: \n{input_desc}".format(
                        input_desc=self.rads_list))

                if self.config['long_term_release']:
                    atm_discharge_C14 = [self.config['annual_discharge_bq_rad_list'][ndx] for ndx, _ in
                                         enumerate(self.rads_list) if _ == 'C-14'][0]
                if self.config['single_plume']:
                    atm_discharge_C14 = [self.config['instantaneous_release_bq_list'][ndx] for ndx, _ in
                                         enumerate(self.rads_list) if _ == 'C-14'][0]

                # for veg
                ingestion_dose_C14 = []
                sum_dose_veg_C14 = 0
                C_air = atm_discharge_C14 * max_dilutfac_for_distance_secperm3
                for veg_type in [self.config['veg_type_list']]:
                    C_pfw = conc_c14_in_terrestrial_plants(C_air, veg_type=veg_type, S_air=0.20)

                    # to convert it in mSv/y
                    dose_veg_C14 = C_pfw * dcfs[ndx] * DID_veg * float(1000)
                    sum_dose_veg_C14 += dose_veg_C14

                # for animal (milk or meat); fish and egg is not yet implemented
                sum_dose_milk_C14 = 0
                sum_dose_meat_C14 = 0
                for animal_prod in self.config['animal_product_list_for_C14']:
                    # get C_tfwt based on animal feed type (one of the input in input.yaml file)

                    C_pfw = conc_c14_in_terrestrial_plants(C_air, veg_type=self.config['animal_feed_type'], S_air=0.20)
                    C_apf = conc_c14_in_terrestrial_animal(C_pfw, animal_feed_type=self.config['animal_feed_type'],
                                                           animal_product=animal_prod, f_c=1)

                    # concentration to dose using dcf for hto and obt
                    if 'DID_milk' in inges_param_dict_adult and any(
                            x in ['cow_milk', 'goat_milk'] for x in self.config['animal_product_list_for_C14']):
                        # to convert it in mSv/y
                        dose_terres_anim_C14 = C_apf * dcfs[ndx] * DID_milk * float(1000)
                        sum_dose_milk_C14 += dose_terres_anim_C14

                    if 'DID_meat' in inges_param_dict_adult and any(
                            x in ['cow_meat', 'goat_meat', 'lamb_meat', 'beef_meat', 'broiler_meat', 'pork_meat'] for x
                            in self.config['animal_product_list_for_C14']):
                        # to convert it in mSv/y
                        dose_terres_anim_C14 = C_apf * dcfs[ndx] * DID_meat * float(1000)
                        sum_dose_meat_C14 += dose_terres_anim_C14

                ingestion_dose_C14.append(
                    ([sum_dose_veg_C14], [sum_dose_milk_C14], [sum_dose_meat_C14]))
        dosepath = 3
        if 'H-3' in self.rads_list and 'C-14' in self.rads_list and len(self.rads_list) > 2:
            # print('np.array(ingestion_dose_tritium):,',np.array(ingestion_dose_tritium).shape, np.array(ingestion_dose_tritium))
            # print('np.array(ingestion_dose_C14):,',np.array(ingestion_dose_C14).shape, np.array(ingestion_dose_C14))
            # print('np.array(ingestion_dose):,',np.array(ingestion_dose).shape, np.array(ingestion_dose))
            # veg, meat, milk ingestion pathway

            ing_dose_stack = np.concatenate(
                (ingestion_dose, np.array(ingestion_dose_tritium).reshape(1, dosepath),
                 np.array(ingestion_dose_C14).reshape(1, dosepath))).T
            self.ingestion_dose_vals = ing_dose_stack

        if 'H-3' in self.rads_list and 'C-14' not in self.rads_list and len(self.rads_list) > 1:
            ing_dose_stack = np.concatenate((ingestion_dose, np.array(ingestion_dose_tritium).reshape(1, dosepath))).T
            self.ingestion_dose_vals = ing_dose_stack

        if 'H-3' not in self.rads_list and 'C-14' in self.rads_list and len(self.rads_list) > 1:
            ing_dose_stack = np.concatenate((ingestion_dose, np.array(ingestion_dose_C14).reshape(1, dosepath))).T
            self.ingestion_dose_vals = ing_dose_stack

        if 'H-3' in self.rads_list and len(self.rads_list) == 1:
            self.ingestion_dose_vals = np.array(ingestion_dose_tritium).reshape(1, dosepath)

        if 'C-14' in self.rads_list and len(self.rads_list) == 1:
            self.ingestion_dose_vals = np.array(ingestion_dose_C14).reshape(1, dosepath)

        if 'H-3' not in self.rads_list and 'C-14' not in self.rads_list and len(self.rads_list) > 0:
            self.ingestion_dose_vals = np.array(ingestion_dose, dtype=object)
        return self.ingestion_dose_vals

    def plumeshine_dose(self, spatial_distance=100):
        """
            Calculates the plume shine dose.

            Input:
                - TJFD_ALL_MISSING_CORR: Triple joint frequency distribution for all stability categories.
                  It has dimension of 6*10*16.
                - rads_list: Radionuclide list (must be provided in array form)
                - X1: Distance from release point to receptor (x direction)
                - HX: Release height (stack)
                - H1: Wind data measured at this height
                - TIM: Sampling time
                - num_days: Number of days used to generate TJFD
                - n: How much sigma required for integral evaluation; default is 1
                - Either long_term_release or single_plume option should be True.
                - Implemented for long_term_release or single_plume.

            Parameters:
                - spatial_distance: Distance from release point to receptor (default=100).

            Returns:
                Plume Shine Dose.
        """
        X1 = spatial_distance
        WSPEED_K = np.array([0.9, 2.4, 4.25, 8.5, 15.5, 24.5, 34.0, 44.5, 56.0, 68.0]) / 3.6  # converted to m/s
        SECWID = 0.39275  # 22.5 degree into radian
        if self.TJFD_ALL_MISSING_CORR is None:
            self.file_preprocessing()
            self.met_data_to_tjfd()
            self.missing_correction()
            self.calm_correction_factor_calc()
            self.dilution_per_sector(spatial_distance)

        # mu=total attenuation coefficient in air (/m); mu_a = energy absorption coefficient in air (/m)
        # NOT PARALLEL
        if not self.config['run_ps_dose_parallel']:
            n_jobs = 1
        # PARALLEL COMPUTATION
        if self.config['run_ps_dose_parallel']:
            n_jobs = -1

        energies, emission_prob, en_dict, neglected_energies = self.gamma_energy_abundaces(
            master_file="library/Dose_ecerman_final.xlsx",
            sheet_name="gamma_energy_radionuclide")
        if len(neglected_energies) > 0:
            logging.getLogger("Related to Plume shine computation:").info(
                "following gamma energies are neglected for plume shine dose computation (in dict format): "
                "{en_dict}".format(
                    en_dict=neglected_energies))

        if len(en_dict) > 0:
            logging.getLogger("Related to Plume shine computation:").info(
                "energy and abundance of radionuclide (in dict format): {en_dict}".format(en_dict=en_dict))
        else:

            print("WARNING: gamma energies for the requested radionuclides are not available in the present "
                  "database library/Dose_ecerman_final.xlsx (sheet: gamma_energy_radionuclide). "
                  "User is requested to explicity add the energies and abundances in "
                  "library/Dose_ecerman_final.xlsx (sheet: gamma_energy_radionuclide) "
                  "file")

        # unit=photons/m2-s
        # integral_stab_cat_energy_wise = []
        # empty list for pure-beta emitter or not-available gamma is filled with [0]
        energies, emission_prob = self.add_zero_energy_for_pure_beta(energies, emission_prob)

        # find k, mu, mu_a and MFP for all the energies of all radionuclides
        k_mu_mua_MFP_dict = self.get_k_mu_mua_MFP(energies)

        # height correction factor for six stab cat
        height_factors = [self.height_correction_factor(sc) for sc in np.arange(1, 7, 1)]

        def get_all_integral_stab_cat_energy_wise_for_all_rad_parallel(X1, energies):
            """
                Calculate the integral stability category energy-wise for all radionuclides in parallel.

                Args:
                    X1: Description of X1.
                    energies (list): List of energies for each radionuclide.

                Returns:
                    list: A list containing integral stability category energy-wise for all radionuclides.
            """
            all_integral_stab_cat_energy_wise = []
            for all_energy_per_rad in energies:
                limit_list_energy_wise = self.get_limit_lists_per_rad_for_all_energies(
                    k_mu_mua_MFP_dict, X1, all_energy_per_rad)
                integral_stab_cat_energy_wise = []
                for edx, energy in enumerate(all_energy_per_rad):
                    limit_list_per_energy = limit_list_energy_wise[edx]
                    # range(6) is range(total_stab_cat)
                    # long term release
                    if self.config['long_term_release']:
                        results = Parallel(n_jobs=n_jobs, verbose=100)(
                            delayed(adgq_sector_average)(ndx, energy, limit_list_per_energy, k_mu_mua_MFP_dict) for
                            ndx, _ in
                            enumerate(range(6)))
                    # single plume release
                    if self.config['single_plume']:
                        results = Parallel(n_jobs=n_jobs, verbose=100)(
                            delayed(adgq_single_plume)(ndx, energy, limit_list_per_energy, k_mu_mua_MFP_dict) for
                            ndx, _ in
                            enumerate(range(6)))
                    integral_stab_cat_energy_wise.append(list(results))
                all_integral_stab_cat_energy_wise.append(integral_stab_cat_energy_wise)
            return all_integral_stab_cat_energy_wise

        if self.config['single_plume']:

            def adgq_single_plume(ndx, energy, limit_list, k_mu_mua_MFP_dict):
                """
                    Calculate the integral using adaptive Gaussian quadrature for a single plume.

                    Args:
                        ndx (int): Index representing stability category.
                        energy: Description of energy.
                        limit_list: Description of limit_list.
                        k_mu_mua_MFP_dict: Description of k_mu_mua_MFP_dict.

                    Returns:
                        float: The calculated integral value.
                """
                stab_cat = ndx + 1
                limit_list = limit_list[ndx]
                X1 = spatial_distance
                k = k_mu_mua_MFP_dict[energy][0]
                mu = k_mu_mua_MFP_dict[energy][1]
                # for single plume
                Y = self.config['Y']
                Z = self.config['Z']
                print('energy, stab_Cat, X1, Y, Z, k, mu:', energy, stab_cat, X1, Y, Z, k, mu, limit_list)

                expo_xyz = lambda x, y, z: (1 / (
                        2 * np.pi * self.sigmay(stab_cat, x)[0] * self.sigmaz(stab_cat, x)[0])) \
                                           * (((1 + (k * mu * np.sqrt((x - X1) ** 2 + (y - Y) ** 2 + (z - Z) ** 2)))
                                               / (4 * np.pi * ((x - X1) ** 2 + (y - Y) ** 2 + (z - Z) ** 2)))) \
                                           * (np.exp(-mu * np.sqrt((x - X1) ** 2 + (y - Y) ** 2 + (z - Z) ** 2))) \
                                           * (np.exp((-0.5 * y ** 2) / (self.sigmay(stab_cat, x)[0] ** 2))) \
                                           * (np.exp(
                    -0.5 * (z - self.release_height) ** 2 / (self.sigmaz(stab_cat, x)[0] ** 2))
                                              + np.exp(
                            -0.5 * (z + self.release_height) ** 2 / (self.sigmaz(stab_cat, x)[0] ** 2)))
                # used adaptive gaussian quadrature
                int_expo_xyz = integrate.tplquad(expo_xyz, limit_list[0], limit_list[1], limit_list[2],
                                                 limit_list[3], limit_list[4], limit_list[5], epsabs=1.49e-05,
                                                 epsrel=1.49e-05)
                return int_expo_xyz[0]

            t = time.time()
            all_integral_stab_cat_energy_wise = get_all_integral_stab_cat_energy_wise_for_all_rad_parallel(X1, energies)
            print('results of parallel plume shine integrals:', all_integral_stab_cat_energy_wise)
            print('time_taken_in_plume_shine_dose_computation:', time.time() - t)

            pl_sh_sectors_list = []
            # for each radionuclide
            for raddx in range(len(energies)):
                pl_sh_sectors = []
                # for each energy value
                for ndx, (energy_rad, abundance_rad) in enumerate(zip(energies[raddx], emission_prob[raddx])):
                    mu_a = k_mu_mua_MFP_dict[energy_rad][2]
                    # looping through six stab cat
                    for i in np.arange(0, 6, 1):
                        plumeshine_dose = 5 * 10 ** (-4) * energy_rad * mu_a * \
                                          all_integral_stab_cat_energy_wise[raddx][ndx][i] * abundance_rad
                        pl_sh_sectors.append(plumeshine_dose)
                pl_sh_sectors = np.array(pl_sh_sectors, dtype=object).reshape(-1, 6).sum(axis=0)
                pl_sh_sectors_list.append(pl_sh_sectors)
            # HERE ALL YEAR CORRESPOND TO ONLY ONE YEAR SO NESTED THE RESULT IN AN ARRAY.
            pl_sh_sectors_list_all_year = np.array(pl_sh_sectors_list)[:, None]

        if self.config['long_term_release']:
            X1 = spatial_distance

            def adgq_sector_average(ndx, energy, limit_list, k_mu_mua_MFP_dict):
                """
                    Calculate the integral using adaptive Gaussian quadrature for sector average.

                    Args:
                        ndx (int): Index representing stability category.
                        energy: Description of energy.
                        limit_list: Description of limit_list.
                        k_mu_mua_MFP_dict: Description of k_mu_mua_MFP_dict.

                    Returns:
                        float: The calculated integral value.
                """
                stab_cat = ndx + 1
                limit_list = limit_list[ndx]
                X1 = spatial_distance
                Y = 0
                Z = self.config['Z']
                k = k_mu_mua_MFP_dict[energy][0]
                mu = k_mu_mua_MFP_dict[energy][1]
                print('energy, stab_Cat, X1, Y, Z, k, mu:', energy, stab_cat, X1, Y, Z, k, mu, limit_list)
                # for sector average
                expo_xyz = lambda x, y, z: (1 / (np.sqrt(2 * np.pi) * x * SECWID * self.sigmaz(stab_cat, x)[0])) \
                                           * ((
                        (1 + (k * mu * np.sqrt((x - X1) ** 2 + (y - Y) ** 2 + (z - Z) ** 2)))
                        / (4 * np.pi * ((x - X1) ** 2 + (y - Y) ** 2 + (z - Z) ** 2)))) \
                                           * (np.exp(
                    -mu * np.sqrt((x - X1) ** 2 + (y - Y) ** 2 + (z - Z) ** 2))) \
                                           * (np.exp(
                    -0.5 * (z - self.release_height) ** 2 / (self.sigmaz(stab_cat, x)[0] ** 2))
                                              + np.exp(
                            -0.5 * (z + self.release_height) ** 2 / (self.sigmaz(stab_cat, x)[0] ** 2)))
                # used adaptive gaussian quadrature
                int_expo_xyz = integrate.tplquad(expo_xyz, limit_list[0], limit_list[1], limit_list[2],
                                                 limit_list[3], limit_list[4], limit_list[5], epsabs=1.49e-03,
                                                 epsrel=1.49e-03)
                return int_expo_xyz[0]

            if not self.config['have_met_data']:
                t = time.time()
                all_integral_stab_cat_energy_wise = get_all_integral_stab_cat_energy_wise_for_all_rad_parallel(X1,
                                                                                                               energies)
                print('results of parallel plume shine integrals:', all_integral_stab_cat_energy_wise)
                print('time_taken_in_plume_shine_dose_computation:', time.time() - t)

                pl_sh_sectors_list = []
                # for each radionuclide
                for raddx in range(len(energies)):
                    pl_sh_sectors = []
                    # for each energy value
                    for ndx, (energy_rad, abundance_rad) in enumerate(zip(energies[raddx], emission_prob[raddx])):
                        mu_a = k_mu_mua_MFP_dict[energy_rad][2]
                        # looping through six stab cat
                        for i in np.arange(0, 6, 1):
                            plumeshine_dose = 5 * 10 ** (-4) * energy_rad * mu_a * \
                                              all_integral_stab_cat_energy_wise[raddx][ndx][i] * abundance_rad
                            pl_sh_sectors.append(plumeshine_dose)
                    pl_sh_sectors = np.array(pl_sh_sectors, dtype=object).reshape(-1, 6).sum(axis=0)
                    pl_sh_sectors_list.append(pl_sh_sectors)
                # HERE ALL YEAR CORRESPOND TO ONLY ONE YEAR SO NESTED THE RESULT IN AN ARRAY.
                pl_sh_sectors_list_all_year = np.array(pl_sh_sectors_list)[:, None]

            # the case where user has provided Met data and long term release
            if self.config['have_met_data']:
                ################## PARALLEL PLUME SHINE CALCULATION ###################################################
                t = time.time()
                all_integral_stab_cat_energy_wise = get_all_integral_stab_cat_energy_wise_for_all_rad_parallel(X1,
                                                                                                               energies)
                print('results of parallel plume shine integrals:', all_integral_stab_cat_energy_wise)
                print('time_taken_in_plume_shine_dose_computation:', time.time() - t)
                ####################### END ############################################################################
                pl_sh_sectors_list_all_year = []
                # for each year
                for each_year in range(len(self.config['excel_sheet_name'])):
                    total_calm_hours = sum(self.TJFD_ALL_MISSING_CORR[each_year][:, 0, :].sum(axis=1))
                    hours_without_calm = (self.num_days[each_year] * self.operation_hours_per_day) - total_calm_hours
                    logging.getLogger("\nplume shine dose calculation").info(
                        "Total Calm Hours: {} hours.".format(total_calm_hours))
                    # for each radionuclide
                    rads_ps_list = []
                    for raddx, (all_energy_per_rad, all_emission_prob) in enumerate(zip(energies, emission_prob)):
                        energywise_ps_list = []
                        for ndx, energy in enumerate(all_energy_per_rad):
                            mu = k_mu_mua_MFP_dict[energy][1]
                            mu_a = k_mu_mua_MFP_dict[energy][2]
                            abundance_rad = all_emission_prob[ndx]
                            ps_factors = 5 * 10 ** (-4) * energy * mu_a * abundance_rad
                            time_unit = 3600 / (hours_without_calm * 60 * 60)
                            ps_dir = []
                            # screen the pure-beta emitter (energy = 0, abundance = 0)
                            if energy > 0 and abundance_rad > 0:
                                # divide each count in TJFD (i.e. a[speed, direc] with speed value for a particular direction
                                abc = np.einsum('ijkl, k->ijkl', [self.TJFD_ALL_MISSING_CORR[each_year][:, 1:, :]],
                                                1 / WSPEED_K[1:])
                                # height correction: height factor shape: (6,)
                                habc = np.einsum('ijkl, j->ijkl', abc, 1 / np.array(height_factors))
                                # sum over speed dimension for particular direction sector
                                ab = np.einsum('ijkl ->ijl', habc)
                                # inte = integral of shape (1,6) for six stab cat
                                inte = np.array(all_integral_stab_cat_energy_wise[raddx][ndx]).reshape(1, 6)
                                # add unit correction and ps_factors
                                inte = inte * time_unit * ps_factors
                                # multiple inte
                                abinte = np.einsum('ijk, ij->ijk', ab, inte)
                                # sum over stability categories to get plumeshine dose across all 16 sectors
                                plumeshine_dose = np.einsum('ijk ->ik', abinte)
                                # store it in a list
                                ps_dir.append(plumeshine_dose)
                            # neglected radionuclide (mostly pure-beta emitter)
                            else:
                                for each in np.zeros(16):
                                    ps_dir.append(each)
                            energywise_ps_list.append(ps_dir)
                        pl_sh_list = np.array(energywise_ps_list).sum(axis=0)
                        rads_ps_list.append(pl_sh_list)
                    pl_sh_sectors_list_all_year.append(rads_ps_list)
        return pl_sh_sectors_list_all_year

    def zeroing_ingestion(self, df, notrelm):
        """
        Zeroes out ingestion values for milk and meat routes for radionuclides where transfer factors are not available.

        Args:
            df (DataFrame): DataFrame containing ingestion dose values for veg/milk/meat route.
            notrelm (list): List of elements for which transfer factor for terrestrial food is not available.

        Returns:
            DataFrame: DataFrame with ingestion values zeroed out for milk and meat routes for radionuclides
                       where transfer factors are not available.
        """
        for ndx, each in enumerate(self.element_list):
            if each in notrelm:
                rad = self.rads_list[ndx]
                # for milk and meat only
                df.loc[rad][1:] = 0
        return df
