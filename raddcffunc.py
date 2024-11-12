import numpy as np
import os
import pandas as pd
from scipy import integrate
import logging
from metfunc import *
from dosefunc import *

os.getcwd()


class RaddcfFunc(DoseFunc):
    """
        Represents a class for generating output.

        Args:
            - device: Description of device.
            - config: Configuration parameters.
            - logdir: Log directory (default is None).

        Attributes:
            - effective_lambda_of_rads: Effective lambda of radionuclides.
            - _device: Description of device.
            - config: Configuration parameters.
            - logdir: Log directory.
            - list_half_life: List of half-lives for radionuclides.
            - speed_dist_array: Speed distribution array.
            - mean_speed_stab_cat_wise: Mean speed stability category wise.
            - lambda_of_rads: Lambda of radionuclides.
            - dilution_factor_sectorwise: Dilution factor sector-wise.
            - dcf_list_submersion_corr: List of dose conversion factors with submersion correction.
            - dcfs_gs_corr: Dose conversion factors with ground shine correction.
            - max_dilution_factor: Maximum dilution factor.
            - ingestion_dose_vals: Ingestion dose values.
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
            - type_rad: Type of radionuclide.

        Raises:
            - ValueError: If required parameters are not provided or if conflicting configurations are set.
    """
    def __init__(self, device, config, log_file_name, logdir=None):
        """
                Initializes an OutputFunc object.

                Args:
                    - device: Description of device.
                    - config: Configuration parameters.
                    - logdir: Log directory (default is None).
        """
        super().__init__(device, config, log_file_name, logdir=None)
        self.effective_lambda_of_rads = None
        self._device = device
        self.config = config
        self.logdir = logdir
        self.list_half_life = None

        self.speed_dist_array = None
        self.mean_speed_stab_cat_wise = None
        self.lambda_of_rads = None
        self.dilution_factor_sectorwise = None
        self.dcf_list_submersion_corr = None
        self.dcfs_gs_corr = None
        if self.config['have_dilution_factor']:
            self.max_dilution_factor = config['list_max_dilution_factor']
        else:
            self.max_dilution_factor = None
        # self.max_dilution_factor = config['max_dilution_factor']
        self.ingestion_dose_vals = None
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

        if self.config['long_term_release'] and self.config['single_plume']:
            raise ValueError("Choose the mode of release: Instantaneous or Continuous (long-term).")

        if not self.config['long_term_release'] and not self.config['single_plume']:
            raise ValueError("Choose the mode of release: Instantaneous or Continuous (long-term).")
        if not self.config['have_dilution_factor']:
            # NEED it no matter you have met data or not to incorporate height correction factor
            self.measurement_height = config['measurement_height']
            if not self.measurement_height:
                raise ValueError("Please provide the measurement height at which meteorological data is measured.")

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

        self.release_height = config['release_height']
        if not self.release_height:
            raise ValueError("Please provide the release height of plume.")

        self.downwind_distances = config['downwind_distances']
        if self.config['plant_boundary'] not in self.downwind_distances:
            self.downwind_distances.append(self.config['plant_boundary'])

        if not self.downwind_distances:
            raise ValueError("Please provide the list of downwind distances (in metre).")

        if self.config['run_dose_computation']:
            self.type_rad = self.config['type_rad']
            self.rads_list = config['rads_list']
            if not self.rads_list:
                raise ValueError("Must provide the list of radionuclides for dose estimation.")

            self.element_list = config['element_list']
            if not self.element_list:
                raise ValueError("Must provide the list of elements for dose estimation.")

            self.age_group = config['age_group']
            if not self.age_group:
                raise ValueError("Please provide the list of age for which dose needs to be calculated.")

            if config['long_term_release']:
                self.annual_discharge_bq_rad_list = config['annual_discharge_bq_rad_list']
                if not self.annual_discharge_bq_rad_list:
                    raise ValueError("Please provide annual discharge info (Bq/year) for each radionuclide.")

    def find_half_life_and_decay_const_radionuclides(self):
        """
            Find the half-life and decay constant for radionuclides.

            This method reads data from the "RadioToxicityMaster.xls" file to find the half-life and decay constant
            (lambda) for each radionuclide in the `rads_list`.

            Returns:
                tuple: A tuple containing two lists: `list_half_life`, which contains the half-life of each radionuclide in
                    seconds, and `lambda_of_rads`, which contains the decay constant (lambda) of each radionuclide in
                    per second.
            Raises:
                ValueError: If the unit of half-life is not recognized. It should be 'a' or 'y' for year, 'd' for days, 'h'
                    for hours, 'm' for minutes, or 's' for seconds.
        """
        # lambda
        xls = pd.ExcelFile("library/RadioToxicityMaster.xls")
        name = pd.read_excel(xls, "Inhalation CED Sv per Bq Public")
        name.dropna(axis=0, how='all', inplace=True)

        # calculate half-life in second
        lambda_of_rads = []
        list_half_life = []
        for rad in self.rads_list:
            search_string = '|'.join([rad])
            df = name[name['Nuclide'].str.contains(search_string, na=False)]
            if str('d') in df['Hal-life'].iloc[0]:
                x = float(df['Hal-life'].iloc[0][:-1])
                x = x * 3600 * 24
                lambda_of_rad = (0.693 / x)
                lambda_of_rads.append(lambda_of_rad)
                # unit second
                list_half_life.append(x)

            elif str('a') in df['Hal-life'].iloc[0] or str('y') in df['Hal-life'].iloc[0]:
                x = float(df['Hal-life'].iloc[0][:-1])
                x = x * 3600 * 24 * 365
                lambda_of_rad = (0.693 / x)
                lambda_of_rads.append(lambda_of_rad)
                list_half_life.append(x)

            elif str('h') in df['Hal-life'].iloc[0]:
                x = float(df['Hal-life'].iloc[0][:-1])
                x = x * 3600
                lambda_of_rad = (0.693 / x)
                lambda_of_rads.append(lambda_of_rad)
                list_half_life.append(x)

            elif str('m') in df['Hal-life'].iloc[0]:
                x = float(df['Hal-life'].iloc[0][:-1])
                x = x * 60
                lambda_of_rad = (0.693 / x)
                lambda_of_rads.append(lambda_of_rad)
                list_half_life.append(x)

            elif str('s') in df['Hal-life'].iloc[0]:
                x = float(df['Hal-life'].iloc[0][:-1])
                lambda_of_rad = (0.693 / x)
                lambda_of_rads.append(lambda_of_rad)
                list_half_life.append(x)

            else:
                raise ValueError('Unit of half-life not recognized. It should be (a or y) for year, d for days, '
                                 'h for hours, m for minutes and s for seconds.')
        # unit: /second
        self.lambda_of_rads = lambda_of_rads
        # unit: second
        self.list_half_life = list_half_life

        return self.list_half_life, self.lambda_of_rads

    def inhalation_dcf_list(self, master_file='library/RadioToxicityMaster.xls',
                            sheet_name='Inhalation CED Sv per Bq Public',
                            age=18):
        """
        Return Dose conversion Factors (Inhalation) specific to age and radionuclide.

        This method reads data from the specified Excel file and sheet to find the dose conversion factors (DCFs) for
        inhalation exposure, considering the age of the public.

        Args:
            master_file (str): The path to the Excel file containing data with DCF.
            sheet_name (str): The name of the sheet in the Excel file containing the respective data.
            age (int): The age of the public.

        Returns:
            numpy.ndarray: An array containing Dose conversion Factors (DCFs) specific to each radionuclide.

        Raises:
            ValueError: If the age of the recipient is not a number.
        """
        xls = pd.ExcelFile(master_file)
        name = pd.read_excel(xls, sheet_name)
        name.dropna(axis=0, how='all', inplace=True)
        type_rad = self.config['type_rad']
        dcfs = []

        for ndx, rad in enumerate(self.rads_list):
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

            df = df[df['Type'].str.contains(search_string_type, na=False)]

            if age > 17:
                dcf = df['e_g_age_g_gt_17a_Sv/Bq'].max()
            elif 12 < age <= 17:
                dcf = df['e_g_age_g_12_17a_Sv/Bq'].max()
            elif 7 < age <= 12:
                dcf = df['e_g_age_g_7_12a_Sv/Bq'].max()
            elif 2 < age <= 7:
                dcf = df['e_g_age_g_2_7a_Sv/Bq'].max()
            elif 1 < age <= 2:
                dcf = df['e_g_age_g_1_2a_Sv/Bq'].max()
            elif age <= 1:
                dcf = df['e_g_age_g_lt_1a_Sv/Bq'].max()
            else:
                raise ValueError('The age of recipient must be a number.')

            dcfs.append(dcf)

        self.inhalation_dcf = np.array(dcfs)
        return self.inhalation_dcf

    def dcf_list_ecerman_ground_shine(self, master_file="library/Dose_ecerman_final.xlsx", sheet_name="surface_dose",
                                      age=18):
        """
        Return Dose conversion Factors (Ground Shine) specific to age and radionuclide.

        This method reads data from the specified Excel file and sheet to find the dose conversion factors (DCFs) for
        ground shine exposure, considering the age of the public.

        Args:
            master_file (str): The path to the Excel file containing data with DCF.
            sheet_name (str): The name of the sheet in the Excel file containing the respective data.
            age (int): The age of the public.

        Returns:
            numpy.ndarray: An array containing Dose conversion Factors (DCFs) specific to each radionuclide.

        Raises:
            ValueError: If the age of the recipient is not a number.
        """
        xls = pd.ExcelFile(master_file)
        name = pd.read_excel(xls, sheet_name)
        name.dropna(axis=0, how='all', inplace=True)
        dcfs = []

        for rad in self.rads_list:
            search_string = '|'.join([rad])
            df = name[name['Nuclide'] == search_string]
            if age > 17:
                dcf = df['Adult'].max()
            elif 12 < age <= 17:
                dcf = df['15-yr-old'].max()
            elif 7 < age <= 12:
                dcf = df['10-yr-old'].max()
            elif 2 < age <= 7:
                dcf = df['5-yr-old'].max()
            elif 1 < age <= 2:
                dcf = df['1-yr-old'].max()
            elif age <= 1:
                dcf = df['Newborn'].max()
            else:
                raise ValueError('The age of recipient must be a number.')

            dcfs.append(dcf)

        self.dcfs_gs = np.array(dcfs)

        return self.dcfs_gs

    def dcf_list_ecerman_ground_shine_include_progeny(self, master_file="library/Dose_ecerman_final.xlsx",
                                                      sheet_name='surface_dose',
                                                      age=18, consider_progeny=True):
        """
            Return Dose conversion Factors (Ground Shine) specific to age and radionuclide, considering progeny contribution.

            This method reads data from the specified Excel file and sheet to find the dose conversion factors (DCFs) for
            ground shine exposure, considering the age of the public. It allows users to choose whether to include the progeny
            contribution in the calculation. DCF values in SRS 19 does not include progeny contribution. This function
            facilitates the computation of DCFs that includes progeny contribution. The details of nuclear decay data
            is obtained from SRS 19 based on ICRP 107.

            Args:
                master_file (str): The path to the Excel file containing data with DCF.
                sheet_name (str): The name of the sheet in the Excel file containing the respective data.
                age (int): The age of the public.
                consider_progeny (bool): A boolean flag indicating whether to include progeny contribution (default is True).

            Returns:
                list of tuples: A list containing tuples of Dose conversion Factors (DCFs) specific to each radionuclide.
                    Each tuple contains two DCFs: one with progeny contribution included and one without.

            Raises:
                ValueError: If the age of the recipient is not a number.
        """
        xls = pd.ExcelFile(master_file)
        name = pd.read_excel(xls, sheet_name)
        name.dropna(axis=0, how='all', inplace=True)

        dcfs = []
        dcf_corr_list = []
        for rad in self.rads_list:

            search_string = '|'.join([rad])
            df = name[name['Nuclide'] == search_string]
            if age > 17:
                dcf_ = df['Adult'].max()
            elif 12 < age <= 17:
                dcf_ = df['15-yr-old'].max()
            elif 7 < age <= 12:
                dcf_ = df['10-yr-old'].max()
            elif 2 < age <= 7:
                dcf_ = df['5-yr-old'].max()
            elif 1 < age <= 2:
                dcf_ = df['1-yr-old'].max()
            elif age <= 1:
                dcf_ = df['Newborn'].max()
            else:
                raise ValueError('The age of recipient must be a number.')
            dcf_corr = dcf_
            daughter_list, frac_yield = self.find_progeny_name_and_yield_f(rad, master_file='library/dcf_corr.xlsx')

            if consider_progeny:
                for d, y in zip(daughter_list, frac_yield):
                    search_string = '|'.join([d])
                    df = name[name['Nuclide'].str.contains(search_string, na=False)]
                    if age > 17:
                        dcf_d = df['Adult'].max()
                    elif 12 < age <= 17:
                        dcf_d = df['15-yr-old'].max()
                    elif 7 < age <= 12:
                        dcf_d = df['10-yr-old'].max()
                    elif 2 < age <= 7:
                        dcf_d = df['5-yr-old'].max()
                    elif 1 < age <= 2:
                        dcf_d = df['1-yr-old'].max()
                    elif age <= 1:
                        dcf_d = df['Newborn'].max()
                    else:
                        raise ValueError('The age of recipient must be a number.')

                    dcf_corr += dcf_d * y
            dcfs.append(dcf_)
            dcf_corr_list.append(dcf_corr)

        if consider_progeny:
            self.dcfs_gs_corr = np.array(dcf_corr_list)
            self.dcfs_gs = np.array(dcfs)

            print('self.dcfs_gs_corr:', self.dcfs_gs_corr, flush=True)
            print('self.dcfs_gs:', self.dcfs_gs, flush=True)
            #print('stacked:', stacked)
            if len(self.dcfs_gs_corr) and len(self.dcfs_gs) > 1:
                stacked = np.vstack((self.dcfs_gs_corr, self.dcfs_gs))
                data_as_tuples = [tuple(inner_list) for inner_list in stacked.T]

                A = stacked[:, 0]
                B = stacked[:, 1]
                return data_as_tuples
            else:
                stacked = [self.dcfs_gs_corr[0], self.dcfs_gs[0]]
                print('stacked', stacked)
                return [stacked]
            #return self.dcfs_gs_corr,  self.dcfs_gs
        else:
            self.dcfs_gs = np.array(dcfs)
            if len(self.dcfs_gs) > 1:

                stacked = np.vstack((self.dcfs_gs, self.dcfs_gs))
                data_as_tuples = [tuple(inner_list) for inner_list in stacked.T]
                A = stacked[:, 0]
                B = stacked[:, 1]
                return data_as_tuples
            else:
                stacked = [self.dcfs_gs[0], self.dcfs_gs[0]]
                return [stacked]

    def reshape_dcfs_with_tritium(self, DCFs_t):
        # Prepare the data to be flattened and structured
        reshaped_dcfs_wd_tritium = []
        index_labels = [1, 18]  # Example index labels
        for i, group in enumerate(DCFs_t):
            for entry in group:
                inhalation = entry[0][0]  # Extract scalar from array
                ground_shine = entry[1][0]  # Access the first element in each list pair
                submersion = entry[2][0]  # Access the first element in each list pair
                ingestion = entry[3].tolist()  # Convert array to list for ingestion values

                # Append a row of flattened data to the data list
                reshaped_dcfs_wd_tritium.append([inhalation, ground_shine, submersion, ingestion])
        return reshaped_dcfs_wd_tritium

    import numpy as np

    def reshape_and_pad_dcf(self, DCFs):
        """
        Reshapes and pads the nested list `DCFs` to a consistent shape [m, 4, n],
        where `n` is the maximum length of any innermost list.

        Parameters:
        DCFs (list): Nested list structure to be reshaped and padded.

        Returns:
        np.array: A NumPy array with shape [m, 4, n].
        """
        # Step 1: Find the maximum length of the innermost lists (n)
        max_n = max(
            max(len(item) if isinstance(item, list) else 1 for item in sublist)
            for outer in DCFs for sublist in outer
        )

        # Step 2: Pad the innermost lists to match the `max_n` dimension
        padded_DCFs = []
        for outer in DCFs:
            padded_outer = []
            for sublist in outer:
                padded_sublist = []
                for item in sublist:
                    # If the item is a list, pad it to `max_n`
                    if isinstance(item, list):
                        padded_item = item + [0] * (max_n - len(item))
                    else:
                        padded_item = [item] + [0] * (max_n - 1)  # If not a list, pad it as a list
                    padded_sublist.append(padded_item)
                padded_outer.append(padded_sublist)
            padded_DCFs.append(padded_outer)

        # Step 3: Convert to a numpy array
        # Ensure all nested lists have homogeneous dimensions for array conversion
        try:
            array_DCFs = np.array(padded_DCFs, dtype=object)
        except ValueError as e:
            print(f"Error: {e}")
            return None

        return array_DCFs

    # Example usage:
    #array_DCFs = reshape_and_pad_dcf(DCFs)

    # Print the final array and its shape
    #print(array_DCFs)
    #print("Shape:", array_DCFs.shape)

    def find_progeny_name_and_yield_f(self, rad, master_file="library/dcf_corr.xlsx"):
        """
            Return the names and fractional yields of progeny radionuclides for a given radionuclide.

            This method reads data from the specified Excel file to find the names and fractional yields of progeny radionuclides
            for a given radionuclide. It allows users to specify a threshold (in seconds) to ignore progeny with half-lives
            shorter than the threshold.

            Args:
                rad (str): The name of the radionuclide.
                master_file (str): The path to the Excel file containing data with progeny radionuclides and fractional yields.

            Returns:
                tuple: A tuple containing two numpy arrays:
                    - An array of strings representing the names of progeny radionuclides.
                    - An array of floats representing the fractional yields of progeny radionuclides.

            Raises:
                ValueError: If the unit of half-life for radionuclides is not recognized.
        """

        xls = pd.ExcelFile(master_file)
        name = pd.read_excel(xls, sheet_name=0,
                             names=['Element', 'Nuclide', 'halflife_and_daughters', 'decaymode_and_yield',
                                    'alpha_energy', 'electron_energy', 'photon_energy', 'total_mev_per_nt'], skiprows=0)
        name.dropna(axis=0, how='all', inplace=True)
        # 30 minute (default)
        if self.config['consider_progeny']:
            ignore_half_life = self.config['ignore_half_life']
        else:
            ignore_half_life = 0
        daughter_list = []
        frac_yield = []
        lambda_of_rads = []
        search_string = '|'.join([rad])
        df = name[name['Nuclide'] == search_string]
        index = df.index[0]
        i = 1
        while name.iloc[[index + i]].isnull().any()['Nuclide']:
            search_string = '|'.join(name.iloc[[index + i]]['halflife_and_daughters'])
            df1 = name[name['Nuclide'].str.contains(search_string, na=False)]
            for _, ndx in zip(df1['Nuclide'].values, df1['Nuclide'].index):
                if str(_) == search_string:
                    hlf = name.iloc[[ndx]]['halflife_and_daughters'].values[0]
                    if str('s') in hlf:
                        x = float(hlf[:-1])

                    elif str('m') in hlf:
                        x = float(hlf[:-1])
                        x = x * 60

                    elif str('d') in hlf:
                        x = float(hlf[:-1])
                        x = x * 3600 * 24

                    elif str('a') in hlf or str('y') in hlf:
                        x = float(hlf[:-1])
                        x = x * 3600 * 24 * 365

                    elif str('h') in hlf:
                        x = float(hlf[:-1])
                        x = x * 3600
                        lambda_of_rad = (0.693 / x)
                        lambda_of_rads.append(lambda_of_rad)

                    else:
                        raise ValueError('the unit of half-life for radionuclides should be in format a (annual),'
                                         'd (day), h (hour), m (minute), s (second)')

                    if x <= ignore_half_life:
                        daughter_list.append(name.iloc[[ndx]]['Nuclide'].values)
                        frac_yield.append(name.iloc[[index + i]]['decaymode_and_yield'].values)
            i += 1

        daughter_list = np.array(daughter_list).flatten()
        frac_yield = np.array(frac_yield).flatten()
        return daughter_list, frac_yield

    def dcf_list_ecerman_submersion(self, master_file='library/Dose_ecerman_final.xlsx', sheet_name='submersion_dose',
                                    age=18):
        """
        Return Dose Conversion Factors (Submersion) specific to age and radionuclide.

        This method retrieves dose conversion factors (DCF) for submersion exposure from the specified Excel file.
        It considers the age of the individual and returns the maximum DCF available for the given radionuclide and age
        group.

        Args:
            master_file (str): The path to the Excel file containing data with DCF for submersion exposure.
            sheet_name (str): The name of the sheet in the Excel file containing the DCF data.
            age (int): The age of the individual for whom DCF is being calculated.

        Returns:
            numpy.ndarray: An array containing the DCF values for submersion exposure specific to each radionuclide
            in the `rads_list`.

        Raises:
            ValueError: If the age of the recipient is not a valid number.
        """
        xls = pd.ExcelFile(master_file)
        name = pd.read_excel(xls, sheet_name)
        name.dropna(axis=0, how='all', inplace=True)
        dcfs = []
        for rad in self.rads_list:
            search_string = '|'.join([rad])
            df = name[name['Nuclide'] == search_string]
            if age > 17:
                dcf = df['Adult'].max()
            elif 12 < age <= 17:
                dcf = df['15-yr-old'].max()
            elif 7 < age <= 12:
                dcf = df['10-yr-old'].max()
            elif 2 < age <= 7:
                dcf = df['5-yr-old'].max()
            elif 1 < age <= 2:
                dcf = df['1-yr-old'].max()
            elif age <= 1:
                dcf = df['Newborn'].max()
            else:
                raise ValueError('The age of recipient must be a number.')

            dcfs.append(dcf)
        self.dcf_list_submersion = np.array(dcfs)

        return self.dcf_list_submersion

    def dcf_list_ecerman_submersion_include_progeny(self, master_file="library/Dose_ecerman_final.xlsx",
                                                    sheet_name="submersion_dose",
                                                    age=18):

        """
        Return Dose Conversion Factors (Submersion) specific to age and radionuclide, including progeny contribution.

        This method computes dose conversion factors (DCF) for submersion exposure from the specified Excel file.
        It accounts for the age of the individual and returns the maximum DCF available for the given radionuclide and age
        group. Additionally, it includes the contribution from progeny radionuclides based on the specified decay data.
        The details of nuclear decay data is obtained from SRS 19 based on ICRP 107.

        Args:
            master_file (str): The path to the Excel file containing data with DCF for submersion exposure.
            sheet_name (str): The name of the sheet in the Excel file containing the DCF data.
            age (int): The age of the individual for whom DCF is being calculated.

        Returns:
            list of tuples: A list containing tuples of DCF values for submersion exposure specific to each radionuclide
            in the `rads_list`, including the contribution from progeny radionuclides if `consider_progeny` is True.
            Each tuple contains two DCF values: the corrected DCF considering progeny and the original DCF.

        Raises:
            ValueError: If the age of the recipient is not a valid number.
        """

        consider_progeny = self.config['consider_progeny']
        xls = pd.ExcelFile(master_file)
        name = pd.read_excel(xls, sheet_name)
        name.dropna(axis=0, how='all', inplace=True)
        dcfs = []
        dcf_corr_list = []
        for rad in self.rads_list:

            search_string = '|'.join([rad])
            df = name[name['Nuclide'] == search_string]
            if age > 17:
                dcf = df['Adult'].max()
            elif 12 < age <= 17:
                dcf = df['15-yr-old'].max()
            elif 7 < age <= 12:
                dcf = df['10-yr-old'].max()
            elif 2 < age <= 7:
                dcf = df['5-yr-old'].max()
            elif 1 < age <= 2:
                dcf = df['1-yr-old'].max()
            elif age <= 1:
                dcf = df['Newborn'].max()
            else:
                raise ValueError('The age of recipient must be a number.')
            daughter_list, frac_yield = self.find_progeny_name_and_yield_f(rad, master_file="library/dcf_corr.xlsx")
            dcf_corr = dcf
            if consider_progeny:
                for d, y in zip(daughter_list, frac_yield):
                    search_string = '|'.join([d])
                    df = name[name['Nuclide'].str.contains(search_string, na=False)]
                    if age > 17:
                        dcf_d = df['Adult'].max()
                    elif 12 < age <= 17:
                        dcf_d = df['15-yr-old'].max()
                    elif 7 < age <= 12:
                        dcf_d = df['10-yr-old'].max()
                    elif 2 < age <= 7:
                        dcf_d = df['5-yr-old'].max()
                    elif 1 < age <= 2:
                        dcf_d = df['1-yr-old'].max()
                    elif age <= 1:
                        dcf_d = df['Newborn'].max()
                    else:
                        raise ValueError('The age of recipient must be a number.')
                    dcf_corr += dcf_d * y

            dcfs.append(dcf)
            dcf_corr_list.append(dcf_corr)

        if consider_progeny:
            self.dcf_list_submersion_corr = np.array(dcf_corr_list)
            print('self.dcf_list_submersion_corr:', self.dcf_list_submersion_corr)

            self.dcf_list_submersion = np.array(dcfs)
            print('self.dcf_list_submersion:', self.dcf_list_submersion)
            if len(self.dcf_list_submersion_corr) > 1 and len(self.dcf_list_submersion) > 1:
                # the output becomes rad specific (corrected, uncorrected DCF)
                stacked = np.vstack((self.dcf_list_submersion_corr, self.dcf_list_submersion))
                data_as_tuples = [tuple(inner_list) for inner_list in stacked.T]
                A = stacked[:, 0]
                B = stacked[:, 1]
                print('A:', A)
                print('B', B)
                return data_as_tuples
            else:
                stacked = [self.dcf_list_submersion_corr[0],
                           self.dcf_list_submersion[0]]
                return [stacked]

        else:
            self.dcf_list_submersion = np.array(dcfs)
            if len(self.dcf_list_submersion) > 1:
                stacked = np.vstack((self.dcf_list_submersion, self.dcf_list_submersion))
                data_as_tuples = [tuple(inner_list) for inner_list in stacked.T]
                A = stacked[:, 0]
                B = stacked[:, 1]
                return data_as_tuples
            else:
                stacked = [self.dcf_list_submersion[0],
                           self.dcf_list_submersion[0]]
                return [stacked]

    def dcf_list_ingestion(self, master_file="library/Dose_ecerman_final.xlsx", sheet_name="ingestion_gsr3", age=18):
        """
        Return Dose Conversion Factors (Ingestion) specific to age and radionuclide.

        This method computes dose conversion factors (DCF) for ingestion exposure from the specified Excel file.
        It accounts for the age of the individual and returns the maximum DCF available for the given radionuclide and age
        group. For tritium (H-3), it returns separate DCF values for HTO and OBT forms.

        Args:
            master_file (str): The path to the Excel file containing data with DCF for ingestion exposure.
            sheet_name (str): The name of the sheet in the Excel file containing the DCF data.
            age (int): The age of the individual for whom DCF is being calculated.

        Returns:
            numpy.ndarray: An array containing DCF values for ingestion exposure specific to each radionuclide in the
            `rads_list`. For H-3, it contains an array with two DCF values for HTO and OBT forms.

        Raises:
            ValueError: If the age of the recipient is not a valid number.
        """
        xls = pd.ExcelFile(master_file)
        name = pd.read_excel(xls, sheet_name)
        name.dropna(axis=0, how='all', inplace=True)
        dcfs = []
        for rad in self.rads_list:
            if rad != 'H-3':
                search_string = '|'.join([rad])
                df = name[name['Nuclide'] == search_string]
                if age > 17:
                    dcf = df['e_g_age_g_gt_17a_Sv/Bq'].max()
                elif 12 < age <= 17:
                    dcf = df['e_g_age_g_12_17a_Sv/Bq'].max()
                elif 7 < age <= 12:
                    dcf = df['e_g_age_g_7_12a_Sv/Bq'].max()
                elif 2 < age <= 7:
                    dcf = df['e_g_age_g_2_7a_Sv/Bq'].max()
                elif 1 < age <= 2:
                    dcf = df['e_g_age_g_1_2a_Sv/Bq'].max()
                elif age <= 1:
                    dcf = df['e_g_age_g_lt_1a_Sv/Bq'].max()
                else:
                    raise ValueError('The age of recipient must be a number.')
                dcfs.append(dcf)
            else:
                dcfs_tritium = []
                for _ in ['HTO', 'OBT']:
                    search_string = '|'.join([_])
                    df = name[name['Nuclide'] == search_string]
                    if age > 17:
                        dcf = df['e_g_age_g_gt_17a_Sv/Bq'].max()
                    elif 12 < age <= 17:
                        dcf = df['e_g_age_g_12_17a_Sv/Bq'].max()
                    elif 7 < age <= 12:
                        dcf = df['e_g_age_g_7_12a_Sv/Bq'].max()
                    elif 2 < age <= 7:
                        dcf = df['e_g_age_g_2_7a_Sv/Bq'].max()
                    elif 1 < age <= 2:
                        dcf = df['e_g_age_g_1_2a_Sv/Bq'].max()
                    elif age <= 1:
                        dcf = df['e_g_age_g_lt_1a_Sv/Bq'].max()
                    else:
                        raise ValueError('The age of recipient must be a number.')
                    dcfs_tritium.append(dcf)
                dcfs.append(dcfs_tritium)

        self.dcfs_ingestion = np.array(dcfs, dtype=object)

        return self.dcfs_ingestion

    def fv_list_ecerman_ingestion(self, master_file="library/Dose_ecerman_final.xlsx", sheet_name="eco_param"):
        """
            Return transfer factors and related parameters for ingestion exposure.

            This method retrieves transfer factors and related parameters for ingestion exposure from the specified Excel file.
            It extracts transfer factors, decay rates, and other parameters for elements relevant to ingestion exposure,
            as listed in SRS 19 Table XI on page 67. For elements not found in the table, it appends them to a list and assigns
            zero values to their parameters.

            Args:
                master_file (str): The path to the Excel file containing data with transfer factors and related parameters
                for ingestion exposure.
                sheet_name (str): The name of the sheet in the Excel file containing the transfer factors and related
                parameters.

            Returns:
                tuple: A tuple containing arrays of transfer factors and related parameters for ingestion exposure. The arrays
                include:
                - lambda_s_per_d_list (numpy.ndarray): Decay rate of radionuclides in soil per day.
                - lambda_w_per_d_list (numpy.ndarray): Decay rate of radionuclides in water per day.
                - f_v1_list (numpy.ndarray): Fraction of radionuclides in soil transferred to the vegetation (roots).
                - f_v2_list (numpy.ndarray): Fraction of radionuclides in soil transferred to the vegetation (leaves).
                - Fm_Milk_d_per_L_list (numpy.ndarray): Transfer factor for radionuclides from feed to milk
                (dose coefficient per liter).
                - Ff_Meat_d_per_kg_list (numpy.ndarray): Transfer factor for radionuclides from feed to meat
                (dose coefficient per kilogram).
                - notransfer_factor_rad (list): List of elements for which transfer factors were not found, and thus,
                their parameters are assigned zero values.

        """

        # from SRS 19 page number: 67, Table XI
        xls = pd.ExcelFile(master_file)
        name = pd.read_excel(xls, sheet_name)
        name.dropna(axis=0, how='all', inplace=True)
        lambda_s_per_d_list = []
        lambda_w_per_d_list = []
        f_v1_list = []
        f_v2_list = []
        Fm_Milk_d_per_L_list = []
        Ff_Meat_d_per_kg_list = []
        notransfer_factor_rad = []
        for rad in self.element_list:
            if rad not in ['H', 'C']:
                search_string = '|'.join([rad])
                df = name[name['Element'].str.contains(search_string, na=False)]
                df = df[df['Element'].str.contains(r'(?:\s|^)%s(?:\s|$)' % search_string)]
                if df.empty:
                    notransfer_factor_rad.append(rad)
                    lambda_s_per_d_list.append(float(0))
                    f_v1_list.append(float(0))
                    f_v2_list.append(float(0))
                    lambda_w_per_d_list.append(float(0))
                    Fm_Milk_d_per_L_list.append(float(0))
                    Ff_Meat_d_per_kg_list.append(float(0))
                else:
                    lspd = df['lambda_s_per_d'].array[0]

                    lambda_s_per_d_list.append(lspd)

                    fv1 = df['Fv1'].array[0]
                    f_v1_list.append(fv1)

                    fv2 = df['Fv2'].array[0]
                    f_v2_list.append(fv2)

                    lwpd = df['lambda_w_per_d'].array[0]
                    lambda_w_per_d_list.append(lwpd)

                    Fm_Milk_d_per_L = df['Fm_Milk_d_per_L'].array[0]
                    Fm_Milk_d_per_L_list.append(Fm_Milk_d_per_L)

                    Ff_Meat_d_per_kg = df['Ff_Meat_d_per_kg'].array[0]
                    Ff_Meat_d_per_kg_list.append(Ff_Meat_d_per_kg)

        return np.array(lambda_s_per_d_list, dtype=float), np.array(lambda_w_per_d_list, dtype=float), np.array(
            f_v1_list, dtype=float), \
               np.array(f_v2_list, dtype=float), np.array(Fm_Milk_d_per_L_list, dtype=float), np.array(
            Ff_Meat_d_per_kg_list, dtype=float), notransfer_factor_rad

    def atten_coeff(self, master_file='library/Dose_ecerman_final.xlsx', sheet_name='mass_attenuation_coeff',
                    energy_rad=None):
        """
            Calculate attenuation coefficients for radiation in air.

            This method calculates attenuation coefficients for radiation in air based on the energy of the radiation.
            It retrieves attenuation coefficient data from the specified Excel file and interpolates the coefficients
            based on the given energy of the radiation. taken from NIST database for air medium:
            https://physics.nist.gov/PhysRefData/XrayMassCoef/ComTab/air.html.

            Args:
                master_file (str): The path to the Excel file containing mass attenuation coefficient data.
                sheet_name (str): The name of the sheet in the Excel file containing the mass attenuation coefficient data.
                energy_rad (float): The energy of the radiation in MeV.

            Returns:
                tuple: A tuple containing the total attenuation coefficient and the energy-dependent attenuation coefficient
                for the given radiation energy. Both coefficients are multiplied by air density in g/cm^3 and 100 to convert
                the unit to /m.

        """
        energy_rad = energy_rad
        xls = pd.ExcelFile(master_file)
        colnames = ['energy', 'total_atten_coeff', 'energy_atten_coeff']
        df = pd.read_excel(xls, sheet_name, names=colnames)
        # energy in MeV; total_atten_coeff and energy_atten_coeff in /m (cm2/g*gm/cm3*100)
        total_atten_coeff = np.interp(energy_rad, df['energy'], df['total_atten_coeff']) * 1.225E-3 * 100
        energy_atten_coeff = np.interp(energy_rad, df['energy'], df['energy_atten_coeff']) * 1.225E-3 * 100
        return total_atten_coeff[0], energy_atten_coeff[0]

    def zyx_lim_for_integral(self, stab_cat, X1, MFP, n=3):
        """
            Calculate spatial limits for triple integral evaluation.

            This method calculates the spatial limits for the triple integral evaluation used in plume dose calculation.
            It considers the stability category, release height, and other parameters to determine the limits for integrating
            over the plume volume.

            Args:
                stab_cat (str): Stability category.
                X1 (float): X-coordinate of the release point.
                MFP (float): Mean free path.
                n (int, optional): Number of standard deviations to consider for the limits. Defaults to 3.

            Returns:
                np.array: An array containing the spatial limits for the triple integral evaluation in the following order:
                    [zin_lim, zf_lim, yin_lim, yf_lim, xin_lim, xf_lim]

        """
        SIGMAZ = self.sigmaz(stab_cat, self.release_height)[0]
        SIGMAY = self.sigmay(stab_cat, self.release_height)[0]
        logging.getLogger("main").info(
            "Sigma Y: {sy}, Sigma Z: {sz}, Stability Category: {sc}".format(sy=SIGMAY, sz=SIGMAZ, sc=stab_cat))

        # find limits for triple integral evaluation
        if self.release_height - (n * SIGMAZ) > 0:
            zin_lim = self.release_height - (n * SIGMAZ)
        else:
            zin_lim = 0
        zf_lim = self.release_height + (n * SIGMAZ)
        yin_lim = -n * SIGMAY
        yf_lim = n * SIGMAY
        '''
        if X1 - (5 * MFP) > 1:
            xin_lim = X1 - (5 * MFP)
        else:
            xin_lim = 1
        '''
        xin_lim = 0
        xf_lim = X1 + (5 * MFP)
        logging.getLogger("Spatial Limits for Plume Dose Calculation").info(
            "zin_lim: {zin_lim}, zf_lim: {zf_lim}, yin_lim: {yin_lim}, yf_lim: {yf_lim}, xin_lim: {xin_lim}, "
            "xf_lim: {xf_lim}" \
                .format(zin_lim=zin_lim, zf_lim=zf_lim, yin_lim=yin_lim, yf_lim=yf_lim, xin_lim=xin_lim, xf_lim=xf_lim))

        return np.array([zin_lim, zf_lim, yin_lim, yf_lim, xin_lim, xf_lim])

    def gamma_energy_abundaces(self, master_file="library/Dose_ecerman_final.xlsx",
                               sheet_name="gamma_energy_radionuclide"):
        """
            Extract gamma energy and abundances data from a specified database.

            This method retrieves gamma energy and emission probability data from a specified database and returns them as lists
            along with a dictionary containing the energy-emission probability pairs. It also identifies and logs neglected
            energies based on certain cutoff criteria. Taken from following database:
            https://www-nds.iaea.org/xgamma_standards/genergies1.htm

            Args:
                master_file (str, optional): Path to the master Excel file containing the data. Defaults to "library/Dose_ecerman_final.xlsx".
                sheet_name (str, optional): Name of the sheet in the Excel file containing the data. Defaults to "gamma_energy_radionuclide".

            Returns:
                tuple: A tuple containing:
                    - energies (list): List of lists containing gamma energies for each radionuclide.
                    - emission_prob (list): List of lists containing emission probabilities for each radionuclide.
                    - en_dict (dict): Dictionary containing energy-emission probability pairs.
                    - neglected_energies (dict): Dictionary containing neglected energies based on cutoff criteria.

        """
        logging.getLogger("main").info("Data source of gamma energy and abundances: {weblk}".format(
            weblk="www-nds.iaea.org/xgamma_standards/genergies1.htm"))
        xls = pd.ExcelFile(master_file)
        colnames = ['nuclide', 'energy_kev', 'std_energy_kev', 'emmission_prob', 'std_emmission_prob', '']
        name = pd.read_excel(xls, sheet_name, header=None, names=colnames)
        name.dropna(axis=0, how='all', inplace=True)
        name = name.iloc[:, :-1]
        energies = []
        emmission_prob = []
        en_dict = {}
        neglected_energies = {}
        for rad in self.rads_list:
            energies_per_rad = []
            emmission_prob_per_rad = []
            search_string = '|'.join([rad])
            df = name[name['nuclide'].str.contains(search_string, na=False)]
            if df.empty:
                emmission_prob_per_rad = [0]
                energies_per_rad = [0]
                print("gamma energy not available for {}. This either means the radionuclide is pure-beta emitter or " 
                      "the data of gamma energy not available in the current database (ref: www-nds.iaea.org/xgamma_standards/genergies1.htm)".format(rad))

                logging.getLogger("gamma energies").info("gamma energy not available for {}. This either means the radionuclide is pure-beta emitter or " 
                      "the data of gamma energy not available in the current database (ref: www-nds.iaea.org/xgamma_standards/genergies1.htm)".format(
                    rad))

            df_e = df['energy_kev'].items()
            df_p = df['emmission_prob'].items()
            for (column_e, content_e), (c, content_p) in zip(df_e, df_p):
                # neglected based on cutoff criteria (below 5 kev or abundance below 1e-03)
                if content_e / 1000 < 0.05 or content_p < 0.001:
                    # converted to MeV unit
                    neglected_energies[rad] = (content_e / 1000, content_p)
                # above 5 kev and abundance above 1e-03
                else:
                    energies_per_rad.append(content_e / 1000)
                    emmission_prob_per_rad.append(content_p)
                    en_dict[content_e / 1000] = content_p

                # converted to MeV unit
                # energies_per_rad.append(content_e / 1000)
                # emmission_prob_per_rad.append(content_p)
                # en_dict[content_e / 1000] = content_p
            energies.append(energies_per_rad)
            emmission_prob.append(emmission_prob_per_rad)
        return energies, emmission_prob, en_dict, neglected_energies

    def deposition_velocity_of_rad(self):
        """
            Calculate the deposition velocity of radionuclides.

            This method calculates the deposition velocity of radionuclides based on various criteria. For certain radionuclides,
            such as tritium (3H), carbon-14 (14C), and non-reactive gases (e.g., helium, neon, argon, krypton, xenon, radon),
            the deposition velocity is assumed to be zero. For some specific reactive gases (e.g., fluorine, chlorine, bromine),
            a fixed deposition velocity of 0.1 m/s is used. For other radionuclides, the deposition velocity is calculated based
            on a specific value, which is consistent with values observed in fallout resulting from certain nuclear accidents.

            Returns:
                list: A list containing the deposition velocity of each radionuclide.

        """
        # deposition velocity; see IAEA SRS 19 (page 27)
        #  total deposition coefficient V_T 1000 m/d be used for screening purposes for deposition of aerosols and reactive gases.
        noble_gas_tritium_c14 = ['H', 'C', 'He', 'Ne', 'Ar', 'Kr', 'Xe', 'Rn']
        noble_gas = ['He', 'Ne', 'Ar', 'Kr', 'Xe', 'Rn']

        list_deposition_vel = []
        for rad in self.element_list:
            # For 3H, 14C and non-reactive gases such as krypton it should be assumed that V_T = 0
            if rad in noble_gas_tritium_c14:
                deposition_vel = 0
                list_deposition_vel.append(deposition_vel)
            elif rad in ['F', 'Cl', 'Br']:
                # unit m/s
                deposition_vel = 0.1
                list_deposition_vel.append(deposition_vel)
            else:
                # earlier 0.01
                # This value of VT is consistent with values for 131I and 137Cs
                # fallout resulting from the accident at the Chernobyl nuclear power station in 1986.
                deposition_vel = float(1000 / 86400)
                list_deposition_vel.append(deposition_vel)

        return list_deposition_vel

    def apply_weathering_correction_gs(self):
        """
            Apply weathering correction to the decay constants for noble gases and selected radionuclides.

            The weathering correction adjusts the decay constants to account for the effects of weathering
            on the radionuclide behavior in the environment. This correction is based on the exposure period
            and specific decay constants for each radionuclide.

            Returns:
            list: A list containing the effective fraction of decay constants after weathering correction
                  for each radionuclide.

            Notes:
            - The weathering correction is applied to the decay constants considering the specific conditions
              for noble gases, iodine (I), technetium (Tc), chlorine (Cl), cesium (Cs), and strontium (Sr).
            - If weathering correction is enabled (config['weathering_corr'] is True), appropriate weathering
              decay constants are added to the original decay constants for the respective radionuclides.
            - The exposure period is converted to seconds before applying the correction.
            - The effective fraction of decay constants after weathering correction is computed using the
              exponential decay formula and the adjusted decay constants.
        """
        # 50 years; plant life # see page 64, SRS 19, the generic value used in the SRS is 30 years.(Table VIII)
        # for decomissioning: 5 years may be choosen
        noble_gas = ['He', 'Ne', 'Ar', 'Kr', 'Xe', 'Rn']
        exposure_period = self.config['exposure_period']
        # convert to second
        exposure_period = exposure_period * 365 * 24 * 3600
        effective_lambda_frac = []
        for idx, lr in enumerate(self.lambda_of_rads):
            rad = self.element_list[idx]
            if self.config['weathering_corr']:
                if rad in noble_gas:
                    weathering_lambda = 0
                    lr_effective = lr + weathering_lambda
                elif rad in ['Tc', 'Cl', 'I']:
                    # see page 66, SRS 19, (Table X); this applies for soil;
                    weathering_lambda = 0.0014 / (24 * 3600)
                    lr_effective = lr + weathering_lambda
                    if rad in ['Tc']:
                        print("For weathering correction for Tc, it is assumed that Tc remains in Tc04- form.")
                    if rad in ['I', 'Cl']:
                        print("For weathering correction for {}, it is assumed that {} remains in anionic form.".format(rad, rad))
                elif rad in ['Cs', 'Sr']:
                    # see page 66, SRS 19, (Table X); this applies for soil;
                    weathering_lambda = 0.00014 / (24 * 3600)
                    lr_effective = lr + weathering_lambda
                else:
                    # see page 66, SRS 19, (Table X); this applies for food crops
                    weathering_lambda = 0
                    lr_effective = lr + weathering_lambda
            else:
                lr_effective = lr
            eff_lambda_per_rad = (1 - np.exp(-lr_effective * exposure_period)) / lr_effective
            effective_lambda_frac.append(eff_lambda_per_rad)
        return effective_lambda_frac

    def ingestion_weathering_correction(self):
        """
            Perform weathering correction for radionuclide concentrations in vegetation.

            This method calculates the weathering correction for radionuclide concentrations in vegetation based on various
            processes including radioactive decay, wash-off by rain or irrigation, surface abrasion, wind action, resuspension,
            tissue aging, leaf fall, herbivore grazing, growth dilution, volatilization, or evaporation. Losses other than
            radioactive decay are described using an aggregated parameter in the form of a first-order rate constant (lw). The
            method applies weathering correction based on the configuration setting and adjusts the effective decay constant
            accordingly.

            Returns:
                list: A list containing the effective decay constant after weathering correction for each radionuclide.
        """
        # PERFORM WEATHERING CORRECTION: The radionuclide concentration in vegetation may be reduced by a variety of
        # processes. These include radioactive decay, wash-off of previously intercepted material by rain or
        # irrigation, surface abrasion and leaf bending from the action of the wind, resuspension, tissue ageing,
        # leaf fall or herbivore grazing, addition of new tissue (growth dilution), volatilization or evaporation.
        # Losses other than radioactive decay are normally described using an aggregated parameter in the form of
        # a first order rate constant lw. Table 7 page 63-64, SRS 19 A default value of lw is given in Table VII for
        # estimating the removal of radionuclides from vegetation based on a half-life of 14 days.

        weathering_corr = self.config['weathering_corr']

        noble_gas_tritium_c14 = ['H', 'C', 'He', 'Ne', 'Ar', 'Kr', 'Xe', 'Rn']
        noble_gas = ['He', 'Ne', 'Ar', 'Kr', 'Xe', 'Rn']
        effective_lambda_of_rads = []

        for idx, lr in enumerate(self.lambda_of_rads):
            rad = self.element_list[idx]
            if weathering_corr:
                if rad in noble_gas_tritium_c14:
                    # For 3H, 14C and non-reactive gases such as krypton it should
                    # be assumed that V_T (dry+wet) = 0
                    lr = lr
                    effective_lambda_of_rads.append(lr)
                elif rad in ['F', 'Cl', 'Br', 'I']:
                    # unit s #check for correctness
                    # Table 7 page 64, SRS 19; UNIT day converted to sec
                    weathering_half_life = 14 * 24 * 3600
                    weathering_lambda = 0.693 / weathering_half_life
                    lr = lr + weathering_lambda
                    effective_lambda_of_rads.append(lr)
                else:
                    # UNIT day converted to sec
                    weathering_half_life = 14 * 24 * 3600
                    weathering_lambda = 0.693 / weathering_half_life
                    lr = lr + weathering_lambda
                    effective_lambda_of_rads.append(lr)
            else:
                lr = lr
                effective_lambda_of_rads.append(lr)

        self.effective_lambda_of_rads = effective_lambda_of_rads

        return self.effective_lambda_of_rads

    def ingestion_weathering_correction_real(self, lambda_w_per_d_list, lambda_s_per_d_list):
        """
            Perform real weathering correction for ingestion using eco-parameters.

            This method calculates the real weathering correction for ingestion using eco-parameters. It adjusts the decay
            constants for radionuclides based on the rate constants for reduction of concentration deposited on plant surfaces
            and decay rates.

            Args:
                lambda_w_per_d_list (list): A list of rate constants for the reduction of concentration deposited on plant
                    surfaces due to processes other than decay (unit: /day).
                lambda_s_per_d_list (list): A list of decay rate constants for radionuclides (unit: /day).

            Returns:
                tuple: A tuple containing two lists:
                    - A list of effective decay constants for ingestion considering weathering correction (unit: /day).
                    - A list of effective decay constants for ingestion without weathering correction (unit: /day).
        """
        # THIS IS WEATHERING CORRECTION FOR INGESTION USING ECO-PARAM ####
        # converted to /day from /second as lambda_w_per_d_list and lambda_s_per_d_list has unit /day
        lambda_i = [each * 24 * 3600 for each in self.lambda_of_rads]
        lambda_eiv_list = []
        lambda_eis_list = []
        for ndx, rad in enumerate(self.element_list):
            if rad not in ['H', 'C']:
                # unit /day; lambda_w_per_d_list: rate constant for reduction of the concentration of material
                # deposited on the plant surface owing to processes apart from decay
                lambda_eiv = lambda_w_per_d_list[ndx] + lambda_i[ndx]
                lambda_eiv_list.append(lambda_eiv)
                # unit /day
                lambda_eis = lambda_s_per_d_list[ndx] + lambda_i[ndx]
                lambda_eis_list.append(lambda_eis)
        return lambda_eiv_list, lambda_eis_list

    def effective_surface_soil_density_rho(self):
        """
            Calculate the effective surface soil density for ingestion dose calculations.

            This method calculates the effective surface soil density based on the soil type specified in the configuration.
            The effective surface soil density is used in ingestion dose calculations.

            Returns:
                tuple: A tuple containing two values:
                    - The effective surface soil density for pasture depth less than 11 cm (unit: g/cm^3).
                    - The effective surface soil density for crop depth greater than or equal to 11 cm (unit: g/cm^3).
        """

        soiltype = self.config['soiltype']
        logging.getLogger("Soil type for ingestion dose calculations").info("Soil: {soil}".format(soil=soiltype))
        if soiltype == 'peatsoil':
            rho_pasture_depth_lt_11 = 50
            rho_crop_depth_ge_11 = 100
        elif soiltype == 'othersoil':
            rho_pasture_depth_lt_11 = 130
            rho_crop_depth_ge_11 = 260
        else:
            raise ValueError('The data for suggested soiltype is not in the database. Option: A) peatsoil B) othersoil')
        return rho_pasture_depth_lt_11, rho_crop_depth_ge_11

    def add_zero_energy_for_pure_beta(self, energies, emission_prob):
        """
            Add zero energy and emission probability for pure beta emitters or cases where gamma energy is not available.

            This method adds a zero energy value and emission probability for pure beta emitters or cases where gamma energy
            is not available in the data. It modifies the energies and emission probabilities lists to include the zero
            energy value.

            Args:
                energies (list): List of gamma energies for radionuclides.
                emission_prob (list): List of corresponding emission probabilities for the gamma energies.

            Returns:
                tuple: A tuple containing two lists:
                    - Modified list of gamma energies, including zero energy where necessary.
                    - Modified list of emission probabilities, including zero values where necessary.
        """
        # unit=photons/m2-s
        # integral_stab_cat_energy_wise = []
        # empty list for pure-beta emitter or not-available gamma is filled with [0]
        energies_mod = []
        emission_prob_mod = []
        for i, j in zip(energies, emission_prob):
            if not i:
                i = [0]
                j = [0]
            emission_prob_mod.append(j)
            energies_mod.append(i)
        return energies_mod, emission_prob_mod

    def zyx_lim_for_integral_single_plume(self, stab_cat, spatial_distance, MFP):
        X1 = float(spatial_distance)
        SIGMAZ = self.sigmaz(stab_cat, X1)[0]
        SIGMAY = self.sigmay(stab_cat, X1)[0]
        logging.getLogger("main").info(
            "Sigma Y: {sigmay}, Sigma Z: {sigmaz}, Stability Category: {sc}".format(sigmay=SIGMAY, sigmaz=SIGMAZ,
                                                                                    sc=stab_cat))

        # find limits for triple integral evaluation
        if self.release_height - (3 * SIGMAZ) > 1:
            zin_lim = self.release_height - (3 * SIGMAZ)
        else:
            zin_lim = 1
        zf_lim = self.release_height + (3 * SIGMAZ)
        yin_lim = -3 * SIGMAY
        yf_lim = 3 * SIGMAY
        # Cite (for MFP): Wang, X. Y., Y. S. Ling, and Z. Q. Shi. "A new finite cloud method for calculating
        # external exposure dose in a nuclear emergency." Nuclear engineering and design 231, no. 2 (2004): 211-216.
        # to solve the problem of singularity the integration is terminated at 1 m downwind distance from
        # the release point
        times = 3
        if X1 - (times * MFP) > 1:
            xin_lim = X1 - (times * MFP)
        else:
            xin_lim = 1
        xf_lim = X1 + (times * MFP)
        logging.getLogger("Spatial Limits for Plume Dose Calculation").info(
            "zin_lim: {zin_lim}, zf_lim: {zf_lim}, yin_lim: {yin_lim}, yf_lim: {yf_lim}, xin_lim: {xin_lim}, "
            "xf_lim: {xf_lim}" \
                .format(zin_lim=zin_lim, zf_lim=zf_lim, yin_lim=yin_lim, yf_lim=yf_lim, xin_lim=xin_lim,
                        xf_lim=xf_lim))

        return np.array([zin_lim, zf_lim, yin_lim, yf_lim, xin_lim, xf_lim])

    def zyx_lim_for_integral_sector_averaged_plume(self, stab_cat, spatial_distance, MFP):
        """
            Calculate the spatial limits for triple integral evaluation for a single plume.

            This method calculates the spatial limits for triple integral evaluation for a single plume based on the stability
            category, spatial distance, and mean free path (MFP).

            Args:
                stab_cat (str): Stability category.
                spatial_distance (float): Spatial distance from the release point.
                MFP (float): Mean free path.

            Returns:
                numpy.ndarray: Array containing the limits for triple integral evaluation in the order [zin_lim, zf_lim,
                yin_lim, yf_lim, xin_lim, xf_lim].
        """
        # self.n = self.config['n']; for y and z direction spreads by n=3 sigma.
        X1 = float(spatial_distance)
        SIGMAZ = self.sigmaz(stab_cat, X1)[0]
        SIGMAY = self.sigmay(stab_cat, X1)[0]
        logging.getLogger("main").info(
            "Sigma Y: {sigmay}, Sigma Z: {sigmaz}, Stability Category: {sc}".format(sigmay=SIGMAY, sigmaz=SIGMAZ,
                                                                                    sc=stab_cat))

        # find limits for triple integral evaluation
        if self.release_height - (3 * SIGMAZ) > 1:
            zin_lim = self.release_height - (3 * SIGMAZ)
        else:
            zin_lim = 1
        zf_lim = self.release_height + (3 * SIGMAZ)
        # next two lines are specific for sector-averaged plume
        yin_lim = -(2 * 3.14 * X1) / (16 * 2)
        yf_lim = (2 * 3.14 * X1) / (16 * 2)
        # Cite (for MFP): Wang, X. Y., Y. S. Ling, and Z. Q. Shi. "A new finite cloud method for calculating
        # external exposure dose in a nuclear emergency." Nuclear engineering and design 231, no. 2 (2004): 211-216.
        # to solve the problem of singularity the integration is terminated at 1 m downwind distance from
        # the release point
        times = 3
        if X1 - (times * MFP) > 1:
            xin_lim = X1 - (times * MFP)
        else:
            xin_lim = 1
        xf_lim = X1 + (times * MFP)
        logging.getLogger("Spatial Limits for Plume Dose Calculation").info(
            "zin_lim: {zin_lim}, zf_lim: {zf_lim}, yin_lim: {yin_lim}, yf_lim: {yf_lim}, xin_lim: {xin_lim}, "
            "xf_lim: {xf_lim}" \
                .format(zin_lim=zin_lim, zf_lim=zf_lim, yin_lim=yin_lim, yf_lim=yf_lim, xin_lim=xin_lim,
                        xf_lim=xf_lim))
        return np.array([zin_lim, zf_lim, yin_lim, yf_lim, xin_lim, xf_lim])

    def atten_coeff(self, master_file="library/Dose_ecerman_final.xlsx", sheet_name="mass_attenuation_coeff",
                    energy_rad=None):
        """
            Calculate the attenuation coefficients for a given energy of radiation in air.

            This function reads attenuation coefficient data from the specified Excel file and sheet, interpolates the values
            for the given energy of radiation, and returns the total attenuation coefficient and energy-dependent attenuation
            coefficient. taken from NIST database for air medium: https://physics.nist.gov/PhysRefData/XrayMassCoef/ComTab/air.html

            Args:
                master_file (str): Path to the Excel file containing the attenuation coefficient data.
                sheet_name (str): Name of the sheet in the Excel file containing the data.
                energy_rad (float): Energy of the radiation in MeV.

            Returns:
                tuple: A tuple containing the total attenuation coefficient and energy-dependent attenuation coefficient
                in /m units.
        """

        energy_rad = energy_rad
        xls = pd.ExcelFile(master_file)
        colnames = ['energy', 'total_atten_coeff', 'energy_atten_coeff']
        df = pd.read_excel(xls, sheet_name, names=colnames)
        # energy in MeV; total_atten_coeff and energy_atten_coeff in /m (cm2/gm * gm/cm3 * 100)
        total_atten_coeff = np.interp(energy_rad, df['energy'], df['total_atten_coeff']) * 1.225E-3 * 100
        energy_atten_coeff = np.interp(energy_rad, df['energy'], df['energy_atten_coeff']) * 1.225E-3 * 100
        return total_atten_coeff[0], energy_atten_coeff[0]

    def get_k_mu_mua_MFP(self, energies):
        """
            Calculate various parameters related to attenuation for a list of energies.

            This function calculates the mass attenuation coefficient (mu), mass energy-absorption coefficient (mu_a),
            mass scattering coefficient (k), and mean free path (MFP) for a list of energies. The values are obtained
            by interpolating data from an Excel file.

            Args:
                energies (list): List of energies in MeV.

            Returns:
                dict: A dictionary containing the calculated parameters for each energy.
        """
        k_mu_mua_MFP_dict = {}
        energies = sum(energies, [])
        for energy in energies:
            mu, mu_a = self.atten_coeff(master_file="library/Dose_ecerman_final.xlsx",
                                        sheet_name="mass_attenuation_coeff",
                                        energy_rad=[energy])
            k = (mu - mu_a) / mu_a
            MFP = 1 / mu
            k_mu_mua_MFP_dict[energy] = (k, mu, mu_a, MFP)
        for key, val in k_mu_mua_MFP_dict.items():
            logging.getLogger("plume_shine").info(
                "For energy {}, values of k={}, mu={}, mua={} and MFP={}'.format(key, val[0], val[1], val[2], val[3])")
            print('For energy {}, values of k={}, mu={}, mua={} and MFP={}'.format(key, val[0], val[1], val[2], val[3]))
        return k_mu_mua_MFP_dict

    def get_limit_lists_per_rad_for_all_energies(self, k_mu_mua_MFP_dict, X1, all_energy_per_rad):
        """
            Calculate limit lists for each energy and stability category combination.

            This function calculates the limit lists for each energy and stability category combination. It iterates over
            all energies and stability categories, obtaining the appropriate parameters from a dictionary, and then calls
            either the zyx_lim_for_integral_sector_averaged_plume or zyx_lim_for_integral_single_plume method to get the
            limit lists.

            Args:
                k_mu_mua_MFP_dict (dict): A dictionary containing parameters calculated for each energy.
                X1 (float): Spatial distance.
                all_energy_per_rad (list): List of energies for each radionuclide.

            Returns:
                numpy.ndarray: A 3D numpy array containing the limit lists for each energy and stability category combination.
                               The shape of the array is (num_energies, 6, 6).
        """
        all_limit_lists_per_rad = []
        for energy in all_energy_per_rad:
            # MFP is obtained using average energy of gamma-emitter; y limit is different in sect-av
            k, mu, mua, MFP = k_mu_mua_MFP_dict[energy]
            for j in np.arange(0, 6, 1):
                stab_cat = j + 1
                # find the integral limits
                if self.config['long_term_release']:
                    limit_list = self.zyx_lim_for_integral_sector_averaged_plume(stab_cat, X1, MFP)
                else:
                    limit_list = self.zyx_lim_for_integral_single_plume(stab_cat, X1, MFP)
                all_limit_lists_per_rad.append(limit_list)
        # make a stack of limit list for six stability categories for each energy. each limit list
        # contains six values
        limit_list_energy_wise = np.array(all_limit_lists_per_rad).reshape(len(all_energy_per_rad), 6, 6)
        return limit_list_energy_wise

    def point_source_dose(self, dist_list=None, activity_curie=1e06, damage_ratio=5e-05,
                          unit='mSv/hr'):
        """
        Calculate the dose rate at specific distances from a point source of radiation.

        This function calculates the dose rate at specific distances from a point source of radiation. It considers the
        gamma energy and yield for each radionuclide, the activity of the source, the distance from the source, and a
        damage ratio for the unanticipated event.

        Args:
            dist_list (list, optional): List of distances in meters. Default is [10, 20, 40, 50, 100, 200, 400, 600, 800, 1000].
            activity_curie (float, optional): Activity of the source in curie. Default is 1e06.
            damage_ratio (float, optional): Damage ratio for the unanticipated event. Default is 5e-05.
            unit (str, optional): Unit of dose rate. Options: 'mSv/hr' or 'mR/hr'. Default is 'mSv/hr'.

        Returns:
            dict: A dictionary containing the dose rate at specific distances for each radionuclide. The keys are the radionuclides,
                  and the values are dictionaries where the keys are the distances and the values are the corresponding dose rates.

        Raises:
            ValueError: If the unit is not recognized. Must be either 'mSv/hr' or 'mR/hr'.
        """

        gamma_energy_list, g_yield_list, en_dict = self.gamma_energy_abundaces()

        if dist_list is None:
            dist_list = [10, 20, 40, 50, 100, 200, 400, 600, 800, 1000]

        if unit == 'mSv/hr':
            dicts = {}
            for ndx, rad in enumerate(self.rads_list()):
                dose_dict = {}
                gamma_energy = gamma_energy_list[ndx]
                g_yield = g_yield_list[ndx]
                for distance in dist_list:
                    dose = 0
                    for ge, y in zip(gamma_energy, g_yield):
                        # unit Sv/hr
                        dose += ((y * 6 * activity_curie * (ge / 1000)) / (distance ** 2)) / (114 * (3.28034 ** 2))
                    # unit mSv/hr
                    dose_after_dr_Corr = dose * damage_ratio * 1000
                    dose_dict[distance] = dose_after_dr_Corr
                dicts.append(dose_dict)

        elif unit == 'mR/hr':
            dicts = {}
            for ndx, rad in enumerate(self.rads_list()):
                dose_dict = {}
                gamma_energy = gamma_energy_list[ndx]
                g_yield = g_yield_list[ndx]
                for distance in dist_list:
                    dose = 0
                    for ge, y in zip(gamma_energy, g_yield):
                        # unit Sv/hr
                        dose += ((y * 6 * activity_curie * (ge / 1000)) / (distance ** 2)) / (3.28034 ** 2)
                    # unit mSv/hr
                    dose_after_dr_Corr = dose * damage_ratio * 1000
                    dose_dict[distance] = dose_after_dr_Corr
                dicts.append(dose_dict)
        else:
            raise ValueError("The unit is not recognized. Must be either 'mSv/hr' or 'mR/hr'. ")

        return dicts

def point_source_dose(gamma_energy=None, g_yield=None, activity_curie=1e06, dist_list=None, damage_ratio=5e-05, unit='mSv/hr'):
    """
        Calculate the dose rate at specific distances from a point source of radiation.

        Args:
            gamma_energy (list, optional): List containing gamma energy (in KeV unit). Default is for Ir-192.
            g_yield (list, optional): Fraction yield/abundance of gamma energy.
            activity_curie (float, optional): Activity of the source in curie. Default is 1e06.
            dist_list (list, optional): List of distances in meters. Default is [10, 20, 40, 50, 100, 200, 400, 600, 800, 1000].
            damage_ratio (float, optional): Damage ratio for the unanticipated event. Default is 5e-05.
            unit (str, optional): Unit of dose rate. Options: 'mSv/hr' or 'mR/hr'. Default is 'mSv/hr'.

        Returns:
            dict: A dictionary containing the dose rate at specific distances.
                  Keys are the distances and values are the corresponding dose rates.

        Raises:
            ValueError: If the unit is not recognized. Must be either 'mSv/hr' or 'mR/hr'.
    """

    if dist_list is None:
        dist_list = [10, 20, 40, 50, 100, 200, 400, 600, 800, 1000]
    if g_yield is None:
        g_yield = [0.0334, 0.2872, 0.2968, 0.8275, 0.4781, 0.03189, 0.04517, 0.082, 0.0534]
    if gamma_energy is None:
        gamma_energy = [205.7943, 295.9565, 308.45507, 316.50618, 468.06885, 484.5751, 588.581, 604.41105,
                        612.46215]
    if unit == 'mSv/hr':
        dose_dict = {}
        for distance in dist_list:
            dose = 0
            for ge, y in zip(gamma_energy, g_yield):
                # unit Sv/hr
                dose += ((y * 6 * activity_curie * (ge / 1000)) / (distance ** 2)) / (114 * (3.28034 ** 2))
            # unit mSv/hr
            dose_after_dr_Corr = dose * damage_ratio * 1000
            dose_dict[distance] = dose_after_dr_Corr
    elif unit == 'mR/hr':
        dose_dict = {}
        for distance in dist_list:
            dose = 0
            for ge, y in zip(gamma_energy, g_yield):
                # unit Sv/hr
                dose += ((y * 6 * activity_curie * (ge / 1000)) / (distance ** 2)) / (3.28034 ** 2)
            # unit mSv/hr
            dose_after_dr_Corr = dose * damage_ratio * 1000
            dose_dict[distance] = dose_after_dr_Corr
    else:
        raise ValueError("The unit is not recognized. Must be either 'mSv/hr' or 'mR/hr'. ")

    return dose_dict






