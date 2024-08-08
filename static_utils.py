import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


class SMetFunc:
    def __init__(self):
        super().__init__(logdir=None)

    @staticmethod
    def file_preprocessing(filepath, sheet_names, start_operation_time=None, end_operation_time=None,
                           column_names=None):
        """
        Process Excel file with meteorological data and returns a 3D array (speed, direction, and stability category).

        Parameters:
        - filepath (str): Path to the Excel file.
        - sheet_names (list): List of sheet names.
        - start_operation_time (str): Start time for operation.
        - end_operation_time (str): End time for operation.
        - column_names (list): List of column names including time, speed, direction, and stability category.

        Returns:
        numpy.ndarray: 3-Dimensional array consisting of meteorological information from the input Excel file.

        """
        # Check if filepath and sheet_names are provided
        if filepath and sheet_names:
            xls = pd.ExcelFile(filepath)
            df_list = []

            # Iterate over each year
            for each_year in range(len(sheet_names)):
                df = pd.read_excel(xls, sheet_names[each_year])

                # Filter data based on operation time
                if start_operation_time and end_operation_time and column_names:
                    if column_names[0] in df.columns:
                        df = df.loc[
                            (df[column_names[0]] >= start_operation_time) & (df[column_names[0]] <= end_operation_time)]

                        # Fill missing values with '999' and reindex
                        f = [c for c in df.columns if c not in column_names[3]]
                        df.loc[:, f] = df.loc[:, f].fillna('999')
                        df = df.reindex(method='nearest')

                        # Handle stability category
                        if isinstance(df[column_names[3]], float) or isinstance(df[column_names[3]], int):
                            df.loc[:, column_names[3]] = df.loc[:, column_names[3]].fillna('9')
                            df[column_names[3]] = [chr(int(x) + 64) for x in df[column_names[3]]]
                        else:
                            df.loc[:, column_names[3]] = df.loc[:, column_names[3]].fillna('I')

                        df = df[column_names[1:]]

                        # Convert stability category to numeric values
                        numeric_stab_cat = []
                        for x in df[column_names[3]]:
                            if isinstance(x, str):
                                numeric_stab_cat.append(ord(x) - 64)
                            else:
                                numeric_stab_cat.append(x)

                        df['stabcat'] = numeric_stab_cat
                        del df[df.columns[-2]]
                        df_list.append(df.to_numpy())

            return np.asarray(df_list, dtype=object)
        else:
            return None

    @staticmethod
    def speed_distribution_list(met_data, quantile=0.90):
        """
        Obtain array of speeds and list of mean speed for all size stability categories.

        Args:
        - met_data (numpy.ndarray): Array of meteorological data.
        - quantile (float): Provide quantile value if the user requires removing data above the requested quantile.
          Default is to keep values within the 90th percentile.

        Returns:
        - numpy.ndarray: Array of speeds for each stability category.
        - list: List of mean speeds for all six stability categories.

        Notes:
        - The shape of met_data should be 1 * num_hours_data * 3 or num_hours_data * 3.

        """

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
        speed_dist_array = np.array(speed_dist_list, dtype=object).reshape(6, 1)

        # Calculate mean speed for each stability category
        mean_speed_stab_cat_wise = [_[0].mean() for _ in speed_dist_array]

        return speed_dist_array, mean_speed_stab_cat_wise

    @staticmethod
    def plot_speed_distribution(speed_dist_array, figname='stability_category_wise_speed_distribution'):
        """
        Plot the speed distribution for each stability category.

        Args:
        - speed_dist_array (numpy.ndarray): Array of speed distribution data for each stability category.
        - figname (str): Name of the output figure file. Default is 'stability_category_wise_speed_distribution'.

        Returns:
        None

        Notes:
        - The speed_dist_array is the array of speed distribution data for six stability categories.
        - The resulting plot is saved as a PNG file with the specified figure name.

        """

        # Create subplots for each stability category
        fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(8, 6))

        # Define legends and colors for each stability category
        legends = ['A', 'B', 'C', 'D', 'E', 'F']
        colors = ['b', 'r', 'g', 'k', 'c', 'y']

        # Iterate over each stability category
        for ndx, (color, ax) in enumerate(zip(colors, axes.flatten())):
            # Compute histogram of speed distribution
            counts, bins = np.histogram(np.array(speed_dist_array[ndx][0]))
            # Plot histogram
            ax.hist(bins[:-1], bins, weights=counts, color=color, alpha=0.7)
            # Set legend, x-axis label, and y-axis label
            ax.legend(legends[ndx], frameon=False)
            ax.set_xlabel('Wind Speed (m/s)')
            ax.set_ylabel('Frequency')

        # Adjust layout and save the plot as a PNG file
        fig.tight_layout()
        plt.savefig('%s.png' % figname, dpi=600)

        return

    @staticmethod
    def met_data_to_tjfd(met_data, sheet_names):
        """
        Process Excel file with meteorological data and returns Triple Joint Frequency Distribution (TJFD).

        Args:
        - met_data (numpy.ndarray): Processed meteorological data from the input Excel file.
        - sheet_names (list): List of sheet names corresponding to each year.

        Returns:
        - met_data: Processed meteorological data from the input Excel file.
        - TJFD_ALL: Triple Joint Frequency Distribution (TJFD) for all years and stability categories.

        Notes:
        - If meteorological data is available, the method processes it and computes the TJFD.
        - The TJFD is calculated separately for each year and stability category based on the meteorological data.
        - The TJFD is computed using histogram2d function with predefined wind speed and direction ranges.
        - The TJFD is adjusted to account for missing data and calm conditions.
        - The method logs information about the preprocessing and TJFD calculation.

        """

        if met_data is not None:
            TJFD_ALL_YEARS = []
            for each_year in range(len(sheet_names)):
                # Count missing data for wind speed, wind direction, and stability category
                counter = 0
                for ndx, each in enumerate(met_data[each_year]):
                    if float(each[1]) > 360 or float(each[2]) > 6:
                        counter += 1

                ISTAB9_missing = np.count_nonzero(met_data[each_year].T[2] == 9)
                WS1_missing = np.count_nonzero(met_data[each_year].T[0] == 999)
                WD1_missing = np.count_nonzero(met_data[each_year].T[1].astype(float) > 360)

                valid_data = len(met_data[each_year]) - counter
                logging.getLogger("report of met data").info(
                    "Year: {year}, speed_data_missing:: {s}, direction_data_missing: {d}, stab_cat_missing: {sc}, "
                    "total invalid data: {inv}, valid data: {valid}" \
                        .format(year=sheet_names[each_year], s=WS1_missing, d=WD1_missing, sc=ISTAB9_missing,
                                inv=counter, valid=valid_data))

                # [0,1.8,3.0,5.5,11.5,19.5,29.5,38.5,50.5,61.5,74.5,87.50]
                WSRANGE = [0, 1.8, 3.0, 5.5, 11.5, 19.5, 29.5, 38.5, 50.5, 61.5, 74.5]

                WDRANGE = [0, 11.25, 33.75, 56.25, 78.75, 101.25, 123.75, 146.25,
                           168.75, 191.25, 213.75, 236.25, 258.75, 281.25, 303.75,
                           326.25, 348.75, 360]

                STABS_ALL = []
                for stab_cat in np.arange(1, 7, 1):
                    STABS = []
                    for ndx, each in enumerate(met_data[each_year]):
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
                    TJFD_ALL.append(Final_TZD)
                TJFD_ALL = np.array(TJFD_ALL)
                logging.getLogger("\ndilution factor calculation").info(
                    "Sheet Name: {year}; TJFD before missing correction: \n{TJFD} ".
                    format(year=sheet_names[each_year], TJFD=TJFD_ALL))
                TJFD_ALL_YEARS.append(TJFD_ALL)

            return met_data, TJFD_ALL_YEARS
        else:
            return None, None
