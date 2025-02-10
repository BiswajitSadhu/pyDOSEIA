import numpy as np
import os
import pandas as pd
from scipy import integrate
import logging
from raddcffunc import *

os.getcwd()


class MetFunc(RaddcfFunc):
    """
    Represents a class for meteorological data processing and dose computation.

    Args:
        - device: Description of device.
        - config: Configuration parameters.
        - logdir: Log directory (default is None).

    Attributes:
        - operation_hours_per_day: Operation hours per day.
        - speed_dist_array: Speed distribution array.
        - mean_speed_stab_cat_wise: Mean speed stability category wise.
        - lambda_of_rads: Lambda of radionuclides.
        - dilution_factor_sectorwise: Dilution factor sector-wise.
        - dcf_list_submersion_corr: List of dose conversion factors with submersion correction.
        - dcfs_gs_corr: Dose conversion factors with ground shine correction.
        - dict_max_dilution_factor: Dictionary of maximum dilution factors.
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
    # def __init__(self):
    def __init__(self, device, config, log_file_name, logdir=None):
        super().__init__(device, config, log_file_name, logdir=None)
        """
            Initializes an MetFunc object.

            Args:
                - device: Description of device.
                - config: Configuration parameters.
                - logdir: Log directory (default is None).
        """
        self.operation_hours_per_day = None
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
            self.dict_max_dilution_factor = config['list_max_dilution_factor']

        else:
            self.dict_max_dilution_factor = None
            self.max_dilution_factor = None

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
        if not self.config['have_dilution_factor']:
            self.measurement_height = config['measurement_height']
            if not self.measurement_height:
                raise ValueError("Please provide the measurement height at which meteorological data is measured.")

        if self.config['long_term_release'] and self.config['single_plume']:
            raise ValueError("Choose the mode of release: Instantaneous or Continuous (long-term).")

        if not self.config['long_term_release'] and not self.config['single_plume']:
            raise ValueError("Choose the mode of release: Instantaneous or Continuous (long-term).")

        self.release_height = config['release_height']
        if not self.release_height:
            raise ValueError("Please provide the release height of plume.")

        self.downwind_distances = config['downwind_distances']
        if self.config['plant_boundary'] not in self.downwind_distances:
            self.downwind_distances.append(self.config['plant_boundary'])

        if not self.downwind_distances:
            raise ValueError("Please provide the list of downwind distances (in metre).")

        if self.config['run_dose_computation']:
            if config['long_term_release']:
                self.annual_discharge_bq_rad_list = config['annual_discharge_bq_rad_list']
                if not self.annual_discharge_bq_rad_list:
                    raise ValueError("Please provide annual discharge info (Bq/year) for each radionuclide.")

            self.rads_list = config['rads_list']
            if not self.rads_list:
                raise ValueError("Must provide the list of radionuclides for dose estimation.")

            self.element_list = config['element_list']
            if not self.element_list:
                raise ValueError("Must provide the list of elements for dose estimation.")

            self.age_group = config['age_group']
            if not self.age_group:
                raise ValueError("Please provide the list of age for which dose needs to be calculated.")

    def file_preprocessing(self):
        """
        Process Excel file with meteorological data and returns a 3D array (speed, direction, and stability category).

        Notes:
        - Excel file with multiple sheets can be passed using the 'sheet' argument within the input.yaml file.
        - It corrects issues with blank entries within Excel, treating them as missing data.
        - The stability category within the Excel sheet can be either 1-6 or A to F.

        Returns:
        numpy.ndarray: 3-Dimensional array consisting of meteorological information from the input Excel file.

        """
        # Check if meteorological data is available
        if self.config['have_met_data']:
            start_operation_time = self.config['start_operation_time']
            end_operation_time = self.config['end_operation_time']
            colnames = self.config['column_names']
            xls = pd.ExcelFile(self.filepath)
            df_list = []

            # Iterate over each year
            for each_year in range(len(self.sheet_names)):
                df = pd.read_excel(xls, self.sheet_names[each_year])

                # Filter data based on operation time
                df = df.loc[(df[colnames[0]] >= start_operation_time) & (df[colnames[0]] <= end_operation_time)]

                # Fill missing values with '999' and reindex
                f = [c for c in df.columns if c not in colnames[3]]
                df.loc[:, f] = df.loc[:, f].fillna(999.0)
                df = df.reindex(method='nearest')

                # Handle stability category
                if isinstance(df[colnames[3]], float) or isinstance(df[colnames[3]], int):
                    df.loc[:, colnames[3]] = df.loc[:, colnames[3]].fillna('9')
                    df.STBCLASS = [chr(int(x) + 64) for x in df.STBCLASS]
                else:
                    df.loc[:, colnames[3]] = df.loc[:, colnames[3]].fillna('I')

                df = df[colnames[1:]]

                # Convert stability category to numeric values
                numeric_stab_cat = []
                for x in df.loc[:, colnames[3]]:
                    if isinstance(x, str):
                        numeric_stab_cat.append(ord(x) - 64)
                    else:
                        numeric_stab_cat.append(x)

                df['stabcat'] = numeric_stab_cat
                del df[df.columns[-2]]
                df_list.append(df.to_numpy())

            self.MET_DATA = np.asarray(df_list, dtype=object)
        else:
            self.MET_DATA = None

        return self.MET_DATA

    def speed_distribution_list(self, quantile=0.90):
        """
        Obtain array of speeds and list of mean speed for all size stability categories.

        Args:
        - quantile (float): Provide quantile value if the user requires removing data above the requested quantile.
          Default is to keep values within the 90th percentile.

        Returns:
        - numpy.ndarray: Array of speeds for each stability category.
        - list: List of mean speeds for all six stability categories.

        Notes:
        - The shape of met_data should be 1 * num_hours_data * 3 or num_hours_data * 3.
        - The met_data can be obtained from the file_preprocessing() function.

        """

        # Perform file preprocessing if MET_DATA is not available
        if self.MET_DATA is None:
            self.file_preprocessing()
            met_data = np.concatenate(self.MET_DATA, axis=0)

            if met_data.ndim == 3:
                x = met_data[0]
            else:
                x = met_data

        else:
            met_data = np.concatenate(self.MET_DATA, axis=0)

            if met_data.ndim == 3:
                x = met_data[0]
            else:
                x = met_data

        # Initialize empty list for storing speed distributions
        speed_dist_list = []

        # Iterate over each stability category
        for stability_cat in np.arange(1, 7, 1):
            speed_dist = x[x[:, 2] == stability_cat].T[0]
            df = pd.Series(speed_dist)
            speed_dist = df[df < df.quantile(quantile)]
            speed_dist_list.append(speed_dist)

        # Convert speed distribution list to numpy array
        self.speed_dist_array = np.array(speed_dist_list, dtype=object).reshape(6, 1)

        # Calculate mean speed for each stability category
        self.mean_speed_stab_cat_wise = [_[0].mean() for _ in self.speed_dist_array]

        return self.speed_dist_array, self.mean_speed_stab_cat_wise

    def plot_speed_distribution(self, figname='stability_category_wise_speed_distribution'):
        """
        Plot the speed distribution for each stability category.

        Args:
        - figname (str): Name of the output figure file. Default is 'stability_category_wise_speed_distribution'.

        Returns:
        None

        Notes:
        - The speed_list is the list of speed distribution data for six stability categories.
        - The file_preprocessing() and speed_distribution_list() methods are called if speed distribution data is not available.
        - The resulting plot is saved as a PNG file with the specified figure name.

        """

        import matplotlib.pyplot as plt

        # Perform file preprocessing and obtain speed distribution data if not available
        if self.speed_dist_array is None:
            self.file_preprocessing()
            self.speed_distribution_list()

        # Create subplots for each stability category
        fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(8, 6))

        # Define legends and colors for each stability category
        legends = ['A', 'B', 'C', 'D', 'E', 'F']
        colors = ['b', 'r', 'g', 'k', 'c', 'y']

        # Iterate over each stability category
        for ndx, (color, ax) in enumerate(zip(colors, axes.flatten())):
            # Compute histogram of speed distribution
            counts, bins = np.histogram(np.array(self.speed_dist_array[ndx][0]))
            # Plot histogram
            ax.hist(bins[:-1], bins, weights=counts, color=color, alpha=0.7)
            # Set legend, x-axis label, and y-axis label
            ax.legend(legends[ndx], frameon=False)
            ax.set_xlabel('Wind Speed (m/s)')
            ax.set_ylabel('Frequency')

        # Adjust layout and save the plot as a PNG file
        fig.tight_layout()
        # Define the file name
        filename = f"{figname}.png"

        # Check if the file exists and remove it
        if os.path.exists(filename):
            os.remove(filename)

        # Save the new figure
        plt.savefig(filename, dpi=600) 
        #plt.savefig('%s.png' % figname, dpi=600)

        return

    def met_data_to_tjfd(self):

        """
        Process Excel file with meteorological data and returns Triple Joint Frequency Distribution (TJFD).

        Returns:
        - MET_DATA: Processed meteorological data from the input Excel file.
        - TJFD_ALL: Triple Joint Frequency Distribution (TJFD) for all years and stability categories.

        Notes:
        - If meteorological data is available, the method processes it and computes the TJFD.
        - The TJFD is calculated separately for each year and stability category based on the meteorological data.
        - The TJFD is computed using histogram2d function with predefined wind speed and direction ranges.
        - The TJFD is adjusted to account for missing data and calm conditions.
        - The method logs information about the preprocessing and TJFD calculation.

        """

        if self.config['have_met_data']:
            if self.MET_DATA is None:
                self.file_preprocessing()

            TJFD_ALL_YEARS = []
            for each_year in range(len(self.sheet_names)):
                # Count missing data for wind speed, wind direction, and stability category
                counter = 0
                for ndx, each in enumerate(self.MET_DATA[each_year]):
                    if float(each[1]) > 360 or float(each[2]) > 6:
                        counter += 1

                ISTAB9_missing = np.count_nonzero(self.MET_DATA[each_year].T[2] == 9)
                WS1_missing = np.count_nonzero(self.MET_DATA[each_year].T[0] == 999)
                WD1_missing = np.count_nonzero(self.MET_DATA[each_year].T[1].astype(float) > 360)

                valid_data = len(self.MET_DATA[each_year]) - counter
                logging.getLogger("report of met data").info(
                    "Year: {year}, speed_data_missing:: {s}, direction_data_missing: {d}, stab_cat_missing: {sc}, "
                    "total invalid data: {inv}, valid data: {valid}" \
                        .format(year=self.sheet_names[each_year], s=WS1_missing, d=WD1_missing, sc=ISTAB9_missing,
                                inv=counter, valid=valid_data))

                # [0,1.8,3.0,5.5,11.5,19.5,29.5,38.5,50.5,61.5,74.5,87.50]
                WSRANGE = [0, 1.8, 3.0, 5.5, 11.5, 19.5, 29.5, 38.5, 50.5, 61.5, 74.5]

                WDRANGE = [0, 11.25, 33.75, 56.25, 78.75, 101.25, 123.75, 146.25,
                           168.75, 191.25, 213.75, 236.25, 258.75, 281.25, 303.75,
                           326.25, 348.75, 360]

                STABS_ALL = []
                for stab_cat in np.arange(1, 7, 1):
                    STABS = []
                    for ndx, each in enumerate(self.MET_DATA[each_year]):
                        if each[2] == stab_cat:
                            STABS.append(each)
                    STABS_ALL.append(STABS)
                TJFD_ALL = []

                # Compute TJFD for each stability category
                for idx, stab_cat_data in enumerate(STABS_ALL):
                    stab_cat_data = np.array(stab_cat_data).astype(float)
                    # WHEN MET DATA DOES NOT HAVE ANY DATA OF ONE STABILITY CATEGORY FOR A YEAR.
                    # for example: ASSUME THE YEAR HAS SEEN ALL STAB CATEGORY EXCEPT F
                    if len(np.array(stab_cat_data).flatten()) == 0:
                        TZD = np.zeros((10, 17))
                    else:
                        WS1 = np.array(stab_cat_data.T[0])
                        WD1 = stab_cat_data.T[1]
                        TZD, WSR, WDR = np.histogram2d(WS1, WD1, bins=[WSRANGE, WDRANGE])
                    Calm = np.array([TZD.T[0] + TZD.T[-1]]).T
                    Final_TZD = np.concatenate((Calm, TZD[:, 1:]), axis=1)[:, :-1]
                    Total_sequence_including_CALM_in_all_directions = Final_TZD.sum(axis=0)
                    # logging.getLogger("Triple Joint Frequency Distribution").info(
                    #    "Total Sequence Including CALM in All Directions:: {tot}".format(
                    #        tot=Total_sequence_including_CALM_in_all_directions))
                    TJFD_ALL.append(Final_TZD)
                TJFD_ALL = np.array(TJFD_ALL)
                logging.getLogger("\ndilution factor calculation").info(
                    "Sheet Name: {year}; TJFD before missing correction: \n{TJFD} ".
                    format(year=self.sheet_names[each_year], TJFD=TJFD_ALL))
                TJFD_ALL_YEARS.append(TJFD_ALL)

            self.TJFD_ALL = TJFD_ALL_YEARS
        # If meteorological data is not available, set TJFD_ALL and MET_DATA to None
        if not self.config['have_met_data']:
            self.TJFD_ALL = None
            self.MET_DATA = None

        return self.MET_DATA, self.TJFD_ALL

    def synthetic_TJFD_for_single_plume(self):
        """
            Generate synthetic Triple Joint Frequency Distribution (TJFD) for a single plume scenario.

            Returns:
            - TJFD_Syn_List: Synthetic TJFD array representing six different conditions for a single plume scenario.

            Notes:
            - This method creates scenarios for an instantaneous release of a plume.
            - Six conditions are tested, each considering one hour of data for a particular stability category.
            - For each condition, only one speed data point is available for the first direction (0 degrees), while all other
              entries for wind speed, direction, and stability categories are set to zero.
            - The generated TJFD_Syn_List represents synthetic TJFDs for the single plume scenario.

        """
        import numpy as np

        # Initialize TJFD_Syn_List with zeros
        TJFD_Syn_List = np.zeros((6, 6, 10, 16))

        # Set one speed data point for the first direction (0 degrees) for each stability category
        for ndx, _ in enumerate(TJFD_Syn_List):
            _[ndx][0][0] = 1
        return TJFD_Syn_List

    def missing_correction(self):
        """
        Perform missing correction for missing or incorrect measured values within meteorological data.

        Method:
        - Calculate the fraction of invalid data by dividing the number of invalid data points by the total valid data points.
        - Multiply each element of the TJFD with this fraction to distribute the loss of invalid data.

        Returns:
        - TJFD_ALL_MISSING_CORR: TJFD with missing correction applied.

        Notes:
        - This method requires meteorological data to be available.
        - It first preprocesses the meteorological data if not already done.
        - If TJFD data is missing, it calculates TJFD using the met_data_to_tjfd method.
        - The missing correction is applied to each TJFD for every year separately.
        - The TJFD_ALL_MISSING_CORR represents TJFD after the missing correction is applied.

        """

        if self.config['have_met_data']:
            if self.MET_DATA is None:
                self.file_preprocessing()
                if self.TJFD_ALL is None:
                    self.met_data_to_tjfd()

        if self.config['have_met_data']:
            TJFD_ALL_MISSING_CORR_LIST = []
            for each_year in range(len(self.sheet_names)):
                counter = 0
                for ndx, each in enumerate(self.MET_DATA[each_year]):
                    if float(each[1]) > 360 or float(each[2]) > 6:
                        counter += 1
                valid_data = len(self.MET_DATA[each_year]) - counter
                TJFD_ALL_MISSING_CORR = []
                for tjfd in self.TJFD_ALL[each_year]:
                    factor = counter / valid_data
                    tjfd = tjfd + (tjfd * factor)
                    TJFD_ALL_MISSING_CORR.append(tjfd.astype(int))
                TJFD_ALL_MISSING_CORR = np.array(TJFD_ALL_MISSING_CORR)
                TJFD_ALL_MISSING_CORR_LIST.append(np.array(TJFD_ALL_MISSING_CORR))
                logging.getLogger("\ndilution factor calculation").info(
                    "Year: {year}, TJFD after missing correction: \n{TJFD_ALL_MISSING_CORR}; \n "
                    "Number of data in six stability categories: \n{total_hours}".
                    format(year=self.sheet_names[each_year], TJFD_ALL_MISSING_CORR=TJFD_ALL_MISSING_CORR,
                           total_hours=[np.sum(each) for each in TJFD_ALL_MISSING_CORR]))
            self.TJFD_ALL_MISSING_CORR = np.array(TJFD_ALL_MISSING_CORR_LIST)
        else:
            self.TJFD_ALL_MISSING_CORR = None
        return self.TJFD_ALL_MISSING_CORR

    def calm_correction_factor_calc(self):
        """
        Calculate calm correction factors using TJFD for calm values (1st speed class) (logic from Hukko Bapat manual).

        Returns:
        - calm_correction_factors: An array of 16 calm correction values for 16 sectors.

        Notes:
        - This method requires meteorological data to be available.
        - If TJFD data with missing correction is not available, it preprocesses the meteorological data and applies missing correction.
        - The calm correction factors are calculated based on the TJFD with missing correction.
        - The calculation follows the logic from the Hukko Bapat manual.
        - The calm correction factors are calculated for each year separately and stored in an array.
        - The calm_correction_factors represent the calm correction values for 16 sectors.

        """

        if self.config['have_met_data']:
            if self.TJFD_ALL_MISSING_CORR is None:
                self.file_preprocessing()
                self.met_data_to_tjfd()
                self.missing_correction()
            calm_correction_factors_list = []
            for each_year in range(len(self.sheet_names)):
                N_0 = sum([np.array(each[0]).sum() for each in np.array(self.TJFD_ALL_MISSING_CORR[each_year])])
                calm_correction_factors = []
                for direc in np.arange(0, 16, 1):
                    N_J = np.array([each[1:] for each in np.array(self.TJFD_ALL_MISSING_CORR[each_year])]).reshape(54,
                                                                                                                   16).sum(
                        axis=0)[
                        direc]
                    N_JL = np.array(
                        [self.TJFD_ALL_MISSING_CORR[each_year][stab_cat][1, direc] for stab_cat in range(6)]).sum()
                    N_L = np.array(self.TJFD_ALL_MISSING_CORR[each_year]).reshape(6, 10, 16).sum(axis=1).sum(axis=0)[1]
                    calm_correction_factor = 1 + ((N_0 * N_JL) / (N_L * N_J))
                    calm_correction_factors.append(calm_correction_factor)
                calm_correction_factors_list.append(calm_correction_factors)

            self.calm_correction_factors = np.array(calm_correction_factors_list)
        else:
            self.calm_correction_factors = None
        return self.calm_correction_factors

    def sigmay(self, stab_cat, X1):
        """
        Compute the standard deviation in the y-direction (sigmay) based on stability category and distance.

        Args:
        - stab_cat (int): Stability category (1-6).
        - X1 (float): Downwind distance (in meters).

        Returns:
        - SIGY (float): Standard deviation in the y-direction.
        - AY (float): Coefficient AY for the stability category.

        Notes:
        - This method computes the standard deviation in the y-direction (sigmay) based on the stability category and the downwind distance.
        - The standard deviation is calculated using the formula: SIGY = AY * (X1 ** 0.9031), where AY is a coefficient specific to the stability category.
        - The stability category ranges from 1 to 6, with 1 being the least stable and 6 being the most stable.
        - The downwind distance X1 should be provided in meters.
        - The sampling time is used to apply a correction factor (CF) to the standard deviation calculation, if applicable.
        - The calculated SIGY and AY values are returned as a tuple.

        """
        # SAMPLING CORR to do
        if self.config['have_met_data']:
            if 60 >= self.sampling_time >= 15:
                CF = (self.sampling_time / 3.0) ** (-0.5)
            elif 240 > self.sampling_time > 60:
                CF = (self.sampling_time / 3.0) ** (-0.4)
        ay_list = [0.3658, 0.2751, 0.2089, 0.1471, 0.1046, 0.0722]
        AY = ay_list[stab_cat - 1]
        SIGY = AY * (X1 ** 0.9031)
        # SIGY = SIGY*CF
        # self.sigma_y = SIGY
        # self.AY = AY
        return SIGY, AY

    def sigmaz(self, stab_cat, X):
        """
            Calculate sigma-Z (SIGZ) for atmospheric dispersion modeling.

            Parameters:
                stab_cat (int): Stability category representing atmospheric stability.
                X (float): Downwind distance parameter.

            Returns:
                tuple: A tuple containing SIGZ (sigma-Z), AZ, Q, and R.
                    - SIGZ (float): Sigma-Z value calculated for the given parameters.
                    - AZ (float): Coefficient for the stability category.
                    - Q (float): Coefficient for the stability category.
                    - R (float): Coefficient for the stability category.

            Raises:
                ValueError: If the downwind distance (X) is not a number.
        """
        X = float(X)
        if X < 100:
            az_list = [0.192, 0.156, 0.116, 0.079, 0.063, 0.053]
            q_list = [0.936, 0.922, 0.905, 0.881, 0.871, 0.814]
            r_list = [0, 0, 0, 0, 0, 0]
            AZ = az_list[stab_cat - 1]
            Q = q_list[stab_cat - 1]
            R = r_list[stab_cat - 1]

        elif 100 <= X <= 1000:
            az_list = [0.00066, 0.038, 0.113, 0.222, 0.211, 0.086]
            q_list = [1.941, 1.149, 0.911, 0.725, 0.678, 0.740]
            r_list = [9.27, 3.3, 0, -1.7, -1.3, -0.35]
            AZ = az_list[stab_cat - 1]
            Q = q_list[stab_cat - 1]
            R = r_list[stab_cat - 1]

        elif X > 1000:
            az_list = [0.00024, 0.055, 0.113, 1.26, 6.73, 18.05]
            q_list = [2.094, 1.098, 0.911, 0.516, 0.305, 0.180]
            r_list = [-9.6, 2.0, 0, -13.0, -34.0, -48.6]
            AZ = az_list[stab_cat - 1]
            Q = q_list[stab_cat - 1]
            R = r_list[stab_cat - 1]

        else:
            raise ValueError('The downwind distance must be a number.')

        SIGZ = (AZ * (X ** Q)) + R
        # self.sigma_z = SIGZ
        return SIGZ, AZ, Q, R

    def master_eq_single_plume(self, SIGMAY, SIGMAZ, factor):
        """
            Calculate the concentration of pollutants using the master equation for a single plume.

            The master equation is taken from Hukko and Bapat's work (eq. 2.5) and calculates the concentration
            of pollutants at a given point based on the input parameters.

            Parameters:
                SIGMAY (float): Standard deviation of the plume in the lateral (crosswind) direction.
                SIGMAZ (float): Standard deviation of the plume in the vertical (height) direction.
                factor (float): Factor to scale the speed of dispersion.

            Returns:
                tuple: A tuple containing pre_expo and expo.
                    - pre_expo (float): Pre-exponential factor in the master equation.
                    - expo (float): Exponential term in the master equation.

            Notes:
                By default, if the configuration is set to maximize concentration along the plume's central line
                on the ground, Y and Z are set to 0. Otherwise, the user can provide custom values for Y and Z.
                Speed is assumed to be unit speed and can be scaled using the factor parameter.

        """

        if self.config['max_conc_plume_central_line_gl']:
            Y = 0
            Z = 0
        else:
            Y = self.config['Y']
            Z = self.config['Z']

        # Assuming unit speed and scaling with factor
        speed = 1
        speed = speed * factor
        pre_expo = 1 / (2 * np.pi * SIGMAY * SIGMAZ * speed)
        expo = (np.exp(-1 * np.divide(np.square(Y), (2 * np.square(SIGMAY))))) * \
               (np.exp(-1 * np.divide(np.square(Z - self.release_height), (2 * np.square(SIGMAZ))))
                + np.exp(-1 * np.divide(np.square(Z + self.release_height), (2 * np.square(SIGMAZ)))))

        return pre_expo, expo

    def master_eq_sector_averaged_plume(self, X1, SIGMAZ):
        """
        Calculate the concentration of pollutants using the master equation for sector-averaged plume.

        The master equation is taken from Hukko and Bapat's work (eq. 2.5) and calculates the concentration
        of pollutants at a given point based on the input parameters.

        Parameters:
            X1 (float): Downwind distance parameter for the sector-averaged plume.
            SIGMAZ (float): Standard deviation of the plume in the vertical (height) direction.

        Returns:
            tuple: A tuple containing pre_expo and expo.
                - pre_expo (float): Pre-exponential factor in the master equation.
                - expo (float): Exponential term in the master equation.

        Notes:
            By default, if the configuration is set to maximize concentration along the plume's central line
            on the ground, Z is set to 0. Otherwise, the user can provide a custom value for Z.
            SECWID represents the angular width of the plume sector in radians.

        """

        SECWID = 0.39275  # 22.5 degree into radian
        if self.config['max_conc_plume_central_line_gl']:
            Z = 0
        else:
            Z = self.config['Z']
        # EQ. 2.32 in Hukko and Bapat (Page 50)
        preexpo = (1 / (np.sqrt(2 * np.pi) * X1 * SECWID * SIGMAZ))
        expo = (np.exp(-1 * np.divide(np.square(Z - self.release_height), (2 * np.square(SIGMAZ))))
                + np.exp(-1 * np.divide(np.square(Z + self.release_height), (2 * np.square(SIGMAZ)))))

        return preexpo, expo

    def dilution_per_sector(self, X1):
        """
            Compute dilution factor per sector using meteorological data.

            Parameters:
                X1 (float): Downwind distance in meters.

            Returns:
                ndarray: Array of dilution factors per sector (for long term release with meteorological data) or per
                stability category (for single plume and long-term plume with without meteorological data).

            Notes:
                The function computes the dilution factor for each sector based on meteorological data and configuration
                settings. It takes into account factors such as stability categories, wind speed, release height, and
                operation hours.

        """

        if self.TJFD_ALL_MISSING_CORR is None and self.config['long_term_release']:
            self.file_preprocessing()
            self.met_data_to_tjfd()
            self.missing_correction()

        # first entry correspond to calm is ignored
        WSPEED_K = np.array([0.9, 2.4, 4.25, 8.5, 15.5, 24.5, 34.0, 44.5, 56.0, 68.0]) / 3.6  # converted to m/s

        if self.config['single_plume']:
            self.TJFD_ALL_MISSING_CORR = self.synthetic_TJFD_for_single_plume()
            logging.getLogger("\ndilution factor calculation").info(
                "Doing Calculation of Dilution Factor for Instantaneous release. TJFD: {TJFD_ALL_MISSING_CORR}".
                format(TJFD_ALL_MISSING_CORR=self.TJFD_ALL_MISSING_CORR))

            kbyq_all_stabcat = []
            for each_stab_cat_tjfd in self.TJFD_ALL_MISSING_CORR:
                kbyq_list = []
                for direc in np.arange(0, 16, 1):
                    SUMKQ_J = 0
                    for i in range(len(each_stab_cat_tjfd)):
                        stab_cat = i + 1
                        SUMNU = 0
                        factor = self.height_correction_factor(stab_cat)
                        for speed_cat, speed in enumerate(WSPEED_K):
                            # for instantaneous release, no need to consider met data, speed is considered 1 m/sec.
                            # later met_data for the area is considered to scale the dilution factor based on mean speed
                            # of respective stability category. no correction for release height.
                            SUMNU += each_stab_cat_tjfd[i][speed_cat, direc]
                        # eq. 2.7 in page 16, hukko and bapat manual
                        SIGMAZ, AZ, Q, R = self.sigmaz(stab_cat, X1)
                        SIGMAY, AY = self.sigmay(stab_cat, X1)
                        # use master equation
                        preexpo, expo = self.master_eq_single_plume(SIGMAY, SIGMAZ, factor)
                        KQIJ = preexpo * expo * SUMNU
                        # assuming 1 data equiv 1 hour is used
                        KQIJ = (KQIJ * 3600) / (60 * 60)
                        SUMKQ_J += KQIJ
                    kbyq_list.append(SUMKQ_J)
                kbyq_all_stabcat.append(kbyq_list)

            # kbyq_all_stabcat = np.array(kbyq_all_stabcat, dtype=object)
            kbyq_sum_over_all_sector = np.array(kbyq_all_stabcat, dtype=object).sum(axis=1)
            self.dilution_factor_sectorwise = kbyq_sum_over_all_sector

            # SCALE THE DILUTION FACTOR (WITH MEAN-SPEED) OR NOT
            # if self.config['scaling_dilution_factor_based_on_met_data_speed_distribution']:
            if self.config['have_met_data']:
                self.speed_dist_array, self.mean_speed_stab_cat_wise = self.speed_distribution_list(quantile=0.90)
                # getting mean speed for six categories and dividing with dilution factor of specific category
                self.dilution_factor_sectorwise = self.dilution_factor_sectorwise / np.array(
                    self.mean_speed_stab_cat_wise)[:, None]
                self.max_dilution_factor = self.dilution_factor_sectorwise.T[0].max()
            else:
                if self.config['like_to_scale_with_mean_speed']:
                    # take first all wind sectors for all the stability category and find maximum
                    # self.max_dilution_factor = self.dilution_factor_sectorwise.T[0].max()
                    self.dilution_factor_sectorwise = self.dilution_factor_sectorwise / np.array(self.config[
                                                                                                     'ask_mean_speed_data'])
                    self.max_dilution_factor = self.dilution_factor_sectorwise.T[0].max()

                else:
                    self.max_dilution_factor = self.dilution_factor_sectorwise.T[0].max()

            # CHECK SHAPE OF DILUTION FACTOR ARRAY
            assert self.dilution_factor_sectorwise.shape == (6, )

        if self.config['long_term_release'] and not self.config['have_dilution_factor'] and not self.config['have_met_data']:
            self.TJFD_ALL_MISSING_CORR = self.synthetic_TJFD_for_single_plume()
            logging.getLogger("\ndilution factor calculation").info(
                "Doing Calculation of Dilution Factor for long-term release under conservative assumptions. TJFD: "
                "\n{TJFD_ALL_MISSING_CORR}".
                format(TJFD_ALL_MISSING_CORR=self.TJFD_ALL_MISSING_CORR))

            kbyq_all_stabcat = []
            # SECWID = 0.39275  # 22.5 degree into radian
            for each_stab_cat_tjfd in self.TJFD_ALL_MISSING_CORR:
                kbyq_list = []
                for direc in np.arange(0, 16, 1):
                    SUMKQ_J = 0
                    for i in range(len(each_stab_cat_tjfd)):
                        stab_cat = i + 1
                        factor = self.height_correction_factor(stab_cat)
                        SUMNU = 0
                        for speed_cat, speed in enumerate(WSPEED_K):
                            # for instantaneous release, no need to consider met data, speed is considered 1 m/sec.
                            # later met_data for the area is considered to scale the dilution factor based on mean speed
                            # of respective stability category. no correction for release height.
                            speed = 1
                            SPEED_IK = speed * factor
                            SUMNU += each_stab_cat_tjfd[i][speed_cat, direc] / SPEED_IK
                        SIGMAZ, AZ, Q, R = self.sigmaz(stab_cat, X1)
                        # uses eq. 2.33 in page 50; Hukko and Bapat manual
                        preexpo, expo = self.master_eq_sector_averaged_plume(X1, SIGMAZ)
                        # preexpo = (2 / (np.sqrt(2 * np.pi) * X1 * SECWID * SIGMAZ))
                        # expo = np.exp((-0.5 * self.release_height * self.release_height) / (SIGMAZ * SIGMAZ))
                        KQIJ = preexpo * expo * SUMNU
                        # assuming 1 data equiv 1 hour is used
                        KQIJ = (KQIJ * 3600) / (60 * 60)
                        SUMKQ_J += KQIJ
                    kbyq_list.append(SUMKQ_J)
                kbyq_all_stabcat.append(kbyq_list)

            kbyq_sum_over_all_sector = np.array(kbyq_all_stabcat, dtype=object).sum(axis=1)
            self.dilution_factor_sectorwise = kbyq_sum_over_all_sector
            # CHECK SHAPE OF DILUTION FACTOR ARRAY
            assert self.dilution_factor_sectorwise.shape == (6, )

            # SCALE THE DILUTION FACTOR (WITH MEAN-SPEED) OR NOT (here no met data available)
            # if self.config['scaling_dilution_factor_based_on_met_data_speed_distribution']:
            if self.config['like_to_scale_with_mean_speed']:
                # take first all wind sectors for all the stability category and find maximum
                # self.max_dilution_factor = self.dilution_factor_sectorwise.T[0].max()
                self.dilution_factor_sectorwise = self.dilution_factor_sectorwise / np.array(
                    self.config['ask_mean_speed_data'])
                self.max_dilution_factor = self.dilution_factor_sectorwise.T[0].max()

            else:
                self.max_dilution_factor = self.dilution_factor_sectorwise.T[0].max()

        if self.config['long_term_release'] and self.config['have_met_data'] and not self.config['have_dilution_factor']:
            start_operation_time = self.config['start_operation_time']
            end_operation_time = self.config['end_operation_time']
            self.operation_hours_per_day = np.absolute(start_operation_time - end_operation_time)
            logging.getLogger("\ndilution factor calculation").info("Doing Calculation of Dilution Factor for "
                                                                    "Continuous release.".format())
            # SECWID = 0.39275  # 22.5 degree into radian
            dilution_factor_sectorwise_all_years = []

            for each_year in range(len(self.sheet_names)):
                total_calm_hours = sum(self.TJFD_ALL_MISSING_CORR[each_year][:, 0, :].sum(axis=1))
                hours_without_calm = (self.num_days[each_year] * self.operation_hours_per_day) - total_calm_hours
                logging.getLogger("\ndilution factor calculation").info("Total Calm Hours: {} hours.".format(total_calm_hours))
                kbyq_list = []
                for direc in np.arange(0, 16, 1):
                    SUMKQ_J = 0
                    for i in range(len(self.TJFD_ALL_MISSING_CORR[each_year])):
                        stab_cat = i + 1
                        SIGMAZ, AZ, Q, R = self.sigmaz(stab_cat, X1)
                        factor = self.height_correction_factor(stab_cat)
                        # uses eq. 2.32 in page 50; Hukko and Bapat manual
                        preexpo, expo = self.master_eq_sector_averaged_plume(X1, SIGMAZ)
                        SUMNU = 0
                        # skip CALM class of speed category
                        for speed_cat, speed in enumerate(WSPEED_K[1:]):
                            SPEED_IK = speed * factor
                            # i = stab_cat
                            # Not considering calm class so, [:, 1:, :]
                            # divide count with corrected speed along a particular sector then sum to get SUMNU
                            SUMNU += self.TJFD_ALL_MISSING_CORR[each_year][:, 1:, :][i][speed_cat, direc] / SPEED_IK
                        KQIJ = preexpo * expo * SUMNU
                        KQIJ = (KQIJ * 3600) / (hours_without_calm * 60 * 60)
                        SUMKQ_J += KQIJ
                    kbyq_list.append(SUMKQ_J)
                kbyq_list = np.array(kbyq_list)

                if self.config['have_met_data'] and self.config['calm_correction']:
                    kbyq_list_calm_corrected = []
                    calm_correction_factors = self.calm_correction_factor_calc()[each_year]
                    for dilution_factor, calm_factor in zip(kbyq_list, calm_correction_factors):
                        kbyq_calm_corrected = dilution_factor * calm_factor
                        kbyq_list_calm_corrected.append(kbyq_calm_corrected)
                    kbyq_list = np.array(kbyq_list_calm_corrected)
                    logging.getLogger("\ndilution factor calculation").info(
                        "'Year = {year} , calm correction factors: {calm_correction_factors}".
                        format(year=self.sheet_names[each_year], calm_correction_factors=calm_correction_factors))

                logging.getLogger("\ndilution factor calculation").info(
                    "Sectorwise Dilution Factors for year: {year};  \n{kbyq_list}".
                    format(year=self.sheet_names[each_year], kbyq_list=kbyq_list))

                dilution_factor_sectorwise_all_years.append(kbyq_list)

            # take average of all wind sectors for all the years and find maximum
            self.dilution_factor_sectorwise = np.array(dilution_factor_sectorwise_all_years).mean(axis=0)
            # CHECK SHAPE OF DILUTION FACTOR ARRAY
            assert self.dilution_factor_sectorwise.shape == (16, )
            if self.config['have_met_data']:
                self.max_dilution_factor = max(self.dilution_factor_sectorwise)
            self.plot_speed_distribution(figname='stability_category_wise_speed_distribution')

        return self.dilution_factor_sectorwise

    def plot_dilution_factor(self, downwind_distance):
        """
            Plot the dilution factor for each wind sector.

            Parameters:
                downwind_distance (float): The downwind distance in meters.

            Returns:
                None

            Notes:
                The function generates a polar bar plot showing the dilution factor for each wind sector.
                It uses the dilution factor data computed previously and visualizes it using matplotlib.
        """

        import matplotlib.font_manager as font_manager
        import matplotlib
        import matplotlib.pyplot as plt

        matplotlib.rcParams.update({'font.size': 10})
        matplotlib.rcParams['font.family'] = ['DejaVu Serif']
        font = font_manager.FontProperties(family='DejaVu Serif',
                                           weight='normal',
                                           style='normal', size=10)

        matplotlib.rc('axes', linewidth=2)
        fig = plt.figure(figsize=(8, 8))

        sectors = ['N', 'NNE', 'NE',
                   'ENE', 'E', 'ESE', 'SE',
                   'SSE', 'S', 'SSW', 'SW', 'WSW',
                   'W', 'WNW', 'NW', 'NNW']

        iN = len(sectors)
        arrCnts = (np.array(self.dilution_factor_sectorwise) / max(self.dilution_factor_sectorwise)) * 100
        theta = np.arange(0, 2 * np.pi, 2 * np.pi / iN)
        width = (2 * np.pi) / iN * 0.9
        bottom = 0
        colors = plt.cm.viridis(self.dilution_factor_sectorwise / max(self.dilution_factor_sectorwise))

        ax = fig.add_axes([0.1, 0.1, 0.75, 0.75], polar=True)
        bars = ax.bar(theta, arrCnts, width=width, color=colors, bottom=bottom, antialiased=True, edgecolor='black')

        rotations = np.rad2deg(theta)

        for x, bar, rotation, label in zip(theta, bars, rotations, sectors):
            lab = ax.text(x, bottom + bar.get_height() + 5, label,
                          ha='left', va='center', rotation=rotation, rotation_mode="anchor", color='green')
        plt.savefig('barplot_dilution_factor_at_%sm' % downwind_distance)

        return

    def test_single_plume_glc_hukkoo(self, dil_fac):
        """
            Test the ground level concentrations computed against Hukkoo and Bapat's results.

            Parameters:
                dil_fac (numpy.ndarray): Array containing the computed ground level concentrations.

            Returns:
                None

            Notes:
                The function compares the computed ground level concentrations with the values provided by Hukkoo and Bapat
                for different downwind distances (100, 200, 300, 500, 700, 1000, 1600, 2000, 4000 meters). It asserts that
                the computed values are close to the provided values, indicating a successful test.
        """

        dil_fac = np.array(dil_fac, dtype=float).T
        hukkoo_glc_100m_single_plume = np.array([[2.2840147833303766e-14, 5.590981589039797e-07, 1.5015465916440028e-05,
                                                  1.8534496171567465e-05, 9.310282356034462e-06, 3.6973162422584492e-06,
                                                  9.083901472613404e-07, 4.64985088543347e-07, 5.809360379016258e-08],
                                                 [5.841151371898279e-22, 1.881938226103609e-09, 8.545950204498517e-07,
                                                  1.2296547461666177e-05, 1.689222515044389e-05, 1.3596757373581253e-05,
                                                  6.948090247268007e-06, 4.718541203801888e-06, 1.2717701082737338e-06],
                                                 [7.956707894464337e-42, 1.091680781377853e-14, 2.634808537371747e-09,
                                                  1.504493947251421e-06, 7.156250389638748e-06, 1.2763487559346726e-05,
                                                  1.1755542375765323e-05, 9.483132536764242e-06,
                                                  3.5384718119021803e-06],
                                                 [1.9691262342948032e-107, 1.7554619140561087e-32,
                                                  2.322138749483794e-18,
                                                  1.645624503785968e-10, 3.9826750086653584e-08, 8.73472028784478e-07,
                                                  4.620161828899931e-06, 6.349041184912812e-06, 6.810849638824183e-06],
                                                 [6.22620696059002e-181, 9.192930646978123e-57, 1.5133024136726587e-31,
                                                  1.0210504603823029e-16, 6.745515294647785e-12, 5.643704797281375e-09,
                                                  4.787367464209608e-07, 1.3403801979441503e-06, 4.725940039623769e-06],
                                                 [0.0, 2.5748671022239412e-139, 1.0642330561676144e-74,
                                                  9.189248563624478e-36, 5.821816484706907e-23, 3.877805752304435e-15,
                                                  5.70587584195314e-10, 8.896730306387878e-09, 5.399181471824446e-07]],
                                                dtype=float).T

        # Dictionary to hold Hukkoo and Bapat's ground level concentrations for different downwind distances
        dict_hukkoo_test_glc = {}
        test_dist = [100, 200, 300, 500, 700, 1000, 1600, 2000, 4000]
        for d, glc in zip(test_dist, hukkoo_glc_100m_single_plume):
            dict_hukkoo_test_glc[d] = glc

        # Dictionary to hold computed ground level concentrations for different downwind distances
        dict_dil_fac = {}
        for d, glc in zip(self.downwind_distances, dil_fac):
            dict_dil_fac[d] = glc

        # Compare computed concentrations with Hukkoo and Bapat's values
        for k in dict_dil_fac.keys():
            if k in test_dist:
                assert np.allclose(dict_dil_fac[k], dict_hukkoo_test_glc[k])
                print("Test passed (Verifying Hukkoo and Bapat results)!!!"
                      " on computing ground level concentration for height 100 m "
                      "with release rate 1 Bq/S, mean speed 1 m/S for downwind distance {} m".format(k))

        return

    def height_correction_factor(self, stab_cat):
        """
            Compute the height correction factor based on the stability category.

            Parameters:
                stab_cat (int): Stability category (1 to 6).

            Returns:
                float: Height correction factor.

            Raises:
                ValueError: If stability category is not an integer from 1 to 6.

            Notes:
                The height correction factor adjusts for the difference in release height and measurement height
                in the atmospheric dispersion calculations. It depends on the stability category of the atmosphere.
                If the release height is less than 10 meters, it's set to 10 meters to avoid dividing by zero or negative values.
                Then, the height correction factor is computed using the adjusted release height and the measurement height.
                If the release height is greater than or equal to 10 meters, the factor is computed directly using the actual release height.

            References:
                - Pasquill, F. (1974). Atmospheric Diffusion (3rd ed.). Horwood.

        """
        if stab_cat < 4:
            AN = 0.2
        elif stab_cat == 4:
            AN = 0.25
        elif stab_cat > 4:
            AN = 0.50
        else:
            raise ValueError('Stability category must be an integer (1 to 6).')

        P = AN / (2.0 - AN)

        if self.release_height < 10:
            release_height = 10
            factor = (release_height / self.measurement_height) ** P
        else:
            factor = (self.release_height / self.measurement_height) ** P
        return factor

    def get_max_dilution_factor(self, dilution_factor_sectorwise_all_distances):
        """
            Get the maximum dilution factor for each downwind distance.

            Parameters:
                dilution_factor_sectorwise_all_distances (list of lists): Dilution factor for all sectors for all distances.

            Returns:
                pandas.DataFrame: DataFrame containing the maximum dilution factor for each downwind distance.

            Notes:
                - If the 'have_dilution_factor' configuration is False, the function calculates the maximum dilution factor for each distance.
                - If 'have_dilution_factor' is True, the function assumes that the dilution factor for all sectors for all distances is provided directly.
                - The returned DataFrame has columns ['Max Dilution Factor'] and index representing the downwind distances.
        """
        # X1 = distance
        # shape: len(dist_list), len(sectors): ex. 11, 16
        # dilution_factor_sectorwise_all_distances is the dilution factor for all sectors for all distances
        if not self.config['have_dilution_factor']:
            m = np.array(dilution_factor_sectorwise_all_distances).max(axis=1)
        else:
            m = np.array(dilution_factor_sectorwise_all_distances)
        max_dil_fac_all_dist = pd.DataFrame(m.T, columns=['Max Dilution Factor'], index=[str(dist) for dist in self.downwind_distances])
        return max_dil_fac_all_dist

    # TO-DO
    def compute_plume_rise_neutral_unstable_cat(self, W0=10, x=100, U=2, D_i=5, D_e=8):
        """
            Compute plume rise under neutral or unstable atmospheric conditions for certain stability categories.

            Parameters:
                W0 (float): Exit velocity (m/sec).
                x (float): Downwind distance (m).
                U (float): Wind speed (m/sec) at stack height.
                D_i (float): Internal stack diameter (m).
                D_e (float): External stack diameter (m).

            Returns:
                float: Plume rise (m).

            Notes:
                This method calculates the plume rise for neutral or unstable atmospheric conditions based on certain assumptions
                and empirical formulas. It is applicable for stability categories A, B, C, and D.

            References:
                - AERB/NF/SG/S-1 Guide; Page 44.
                - IAEA-TECDOC-379.
        """
        C = 3 * (1.5 - (W0 / U)) * D_e
        # Dh = Plume rise
        Dh_a = 1.44 * D_i * ((W0 / U) ** (2 / 3)) * ((x / D_i) ** (1 / 3)) - C

        Dh_b = 3 * D_i * (W0 / U)

        Dh = min([Dh_a, Dh_b])
        return Dh

    # TO-DO
    def compute_plume_rise_stable_cat(self, W0=10, U=2, D_i=5):
        """
            Compute plume rise under stable atmospheric conditions for certain stability categories.

            Parameters:
                W0 (float): Exit velocity (m/sec).
                U (float): Wind speed (m/sec) at stack height.
                D_i (float): Internal stack diameter (m).

            Returns:
                float: Plume rise (m).

            Notes:
                This method calculates the plume rise for stable atmospheric conditions based on certain assumptions and
                empirical formulas. It is applicable for stability categories E and F.

            References:
                - AERB/NF/SG/S-1 Guide; Page 44.
                - IAEA-TECDOC-379.
        """
        # Fm is the momentum flux parameter
        Fm = (W0 ** 2) * ((D_i / 2) ** 2)
        # For E stability class,
        S = 8.7E-4
        # For F stability class
        S = 1.75E-3

        Dh = 4 * ((Fm / S) ** (1 / 4))
        Dh = 1.5 * ((S) ** (-1 / 6)) * ((Fm / U) ** (1 / 3))
        return Dh

    # TO-DO
    def building_wake_effect_gifford(self, chi_by_q_original, A, U, Q=1):
        """
            Calculate the wake effect on the plume dispersion behind a building using Gifford's method.

            Parameters:
                chi_by_q_original (float): Original value of chi/Q, where chi is the plume centerline concentration and Q is the emission rate.
                A (float): Cross-sectional area of the building.
                U (float): Mean wind speed.
                Q (float, optional): Emission rate. Default is 1.

            Returns:
                float: Corrected value of chi/Q accounting for building wake effect.

            References:
                - AERB/NF/SG/S-1 Guide; Page 44
                - IAEA-TECDOC-379

            Notes:
                - This function calculates the wake effect on the plume dispersion behind a building based on Gifford's method.
                - The wake effect is determined by the cross-sectional area of the building (A), mean wind speed (U), and emission rate (Q).
                - The calculated value of chi/Q is corrected to ensure it does not fall below one-third of the original value.

        """
        # A = cross-sectional area of the building
        C = 0.5
        # U = mean_speed
        chi_comp = Q / ((C * A) + np.pi * self.sigmay * self.sigmaz) * U
        chi_by_q_comp = chi_comp / Q
        if chi_by_q_comp < (1 / 3) * chi_by_q_original:
            chi_by_q_comp = (1 / 3) * chi_by_q_original

        return chi_by_q_comp
