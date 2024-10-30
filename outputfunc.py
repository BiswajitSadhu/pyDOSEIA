import numpy as np
import os
import pandas as pd
from scipy import integrate
import logging
from metfunc import *
from raddcffunc import *
from dosefunc import *
import csv
import itertools
import _pickle as cPickle
os.getcwd()


# metfunc = MetFunc(config=config, device=args.device, logdir=args.logdir)

class OutputFunc(MetFunc):
    """
        Represents a class for generating output.

        Args:
            - device: Description of device.
            - config: Configuration parameters.
            - logdir: Log directory (default is None).

        Attributes:
            - max_dil_fac_all_dist: Maximum dilution factor for all distances.
            - operation_hours_per_day: Operation hours per day.
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

        Raises:
            - ValueError: If required parameters are not provided or if conflicting configurations are set.
    """
    def __init__(self, device, config, logdir=None):
        super().__init__(device, config, logdir=None)
        self.max_dil_fac_all_dist = None
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

        logging.getLogger("main").info("input description:: \n{input_desc}".format(input_desc=self.config))

    def output_to_txt(self, dilution_factor_sectorwise_all_distances, filename, DCFs=None,
                      DOSES=None, INGESTION_DOSES=None, PLUME_DOSES=None):
        """
            Save results to a text file.

            Parameters:
                dilution_factor_sectorwise_all_distances (list): List of dilution factors for different wind sectors and distances.
                filename (str): Name of the output text file.
                DCFs (list, optional): Dose conversion factors. Defaults to None.
                DOSES (list, optional): Doses computed for different scenarios. Defaults to None.
                INGESTION_DOSES (list, optional): Ingestion doses. Defaults to None.
                PLUME_DOSES (list, optional): Plume shine doses. Defaults to None.

            Returns:
                None

            Notes:
                This method writes the results of the dose computation and other relevant information to a text file.

        """
        # SAVING RESULTS IN AN OUTPUT FILE
        f = open(filename, 'w')
        # unit: metre
        plant_boundary = self.config['plant_boundary']
        pboundary = [str(plant_boundary) + 'm']

        letters = {
            'A': [' *** ', '*   *', '*****', '*   *', '*   *'],
            'B': ['**** ', '*   *', '**** ', '*   *', '**** '],
            'C': [' ****', '*    ', '*    ', '*    ', ' ****'],
            'D': ['**** ', '*   *', '*   *', '*   *', '**** '],
            'E': ['*****', '*    ', '*****', '*    ', '*****'],
            'F': ['*****', '*    ', '***  ', '*    ', '*    '],
            'G': [' ****', '*    ', '*  **', '*   *', ' ****'],
            'H': ['*   *', '*   *', '*****', '*   *', '*   *'],
            'I': ['*****', '  *  ', '  *  ', '  *  ', '*****'],
            'J': ['*****', '    *', '    *', '*   *', ' *** '],
            'K': ['*   *', '*  * ', '**   ', '*  * ', '*   *'],
            'L': ['*    ', '*    ', '*    ', '*    ', '*****'],
            'M': ['*   *', '** **', '* * *', '*   *', '*   *'],
            'N': ['*   *', '**  *', '* * *', '*  **', '*   *'],
            'O': [' *** ', '*   *', '*   *', '*   *', ' *** '],
            'P': ['**** ', '*   *', '**** ', '*    ', '*    '],
            'Q': [' *** ', '*   *', '*   *', '*  **', ' ** *'],
            'R': ['**** ', '*   *', '**** ', '*  * ', '*   *'],
            'S': [' ****', '*    ', '**** ', '    *', '**** '],
            'T': ['*****', '  *  ', '  *  ', '  *  ', '  *  '],
            'U': ['*   *', '*   *', '*   *', '*   *', ' *** '],
            'V': ['*   *', '*   *', '*   *', ' * * ', '  *  '],
            'W': ['*   *', '*   *', '* * *', '** **', '*   *'],
            'X': ['*   *', ' * * ', '  *  ', ' * * ', '*   *'],
            'Y': ['*   *', ' * * ', '  *  ', '  *  ', '  *  '],
            'Z': ['*****', '   * ', '  *  ', ' *   ', '*****'],
        }

        string = "pyDOSEIA"

        f.write('\n')
        f.write(
            '###############################################################################################################################\n')
        # f.write('\n')
        for i in range(5):
            for word in range(len(string)):
                current_word = string[word].upper()
                if word == len(string) - 1:
                    print(letters[current_word][i], file=f)

                else:
                    print(letters[current_word][i], end='  ', file=f)

        pd.set_option('display.float_format', '{:.3e}'.format)

        ###CODE INTRO##################
        f.write(
            '###############################################################################################################################')
        f.write('#\n#\n#')
        f.write(
            '#pyDOSEIA: python-based software package for radiation dose computation in case of atmospheric release                                    ')
        f.write('#\n#\n#')
        f.write(
            '#Code written by: Dr. Biswajit Sadhu, SSS, HPD, HS & E Group (email: bsadhu@barc.gov.in)                                      ')
        f.write('#\n#\n#')
        f.write(
            '#Data compilation and benchmarking: Shri Tanmay Sarkar and Dr. Biswajit Sadhu                                                 ')
        f.write('#\n#\n#')
        f.write(
            '#Conceptualization and supervision: Dr. Biswajit Sadhu, Shri Tanmay Sarkar, Dr. S. Anand, Shri KapilDeo Singh, Dr. M. S. Kulkarni')
        f.write('#\n#\n#')
        f.write(
            '###############################################################################################################################')
        f.write('\n\n')
        f.write('########## YOUR INPUT ##############')
        f.write('\n')
        # pd.options.display.max_colwidth = 90
        f.write("\n".join("{}:\t{}".format(k, v) for k, v in self.config.items()))
        f.write(
            '\nN.B. The height and distances are in metre unit. If given, the unit of ignore_half_life and exposure_period are second and year, respectively.')
        # df_inp = pd.DataFrame(self.config.items(), columns=['input', 'value'])
        # df_inp.style.set_properties(**{'text-align': 'left'})
        # f.write(df_inp.__repr__())
        # f.write(str(self.config))
        f.write('\n')



        sectors = ['N', 'NNE', 'NE',
                   'ENE', 'E', 'ESE', 'SE',
                   'SSE', 'S', 'SSW', 'SW', 'WSW',
                   'W', 'WNW', 'NW', 'NNW']
        stab_cat_name = np.array(['A', 'B', 'C', 'D', 'E', 'F'])
        ###PRINT Dilution factors##################
        if self.config['single_plume']:

            if not self.config['have_dilution_factor']:
                # .sum(axis=2)
                df_dilution_factor = pd.DataFrame(
                    np.array(dilution_factor_sectorwise_all_distances, dtype=object).T,
                    columns=[str(dist) + 'm' for dist in self.downwind_distances],
                    index=stab_cat_name)

                # TESTING THE RESULT OF GROUND LEVEL CONCENTRATION WITH HUKKOO (PAGE 98)
                # CONDITIONl Q=1 Bq/Sec, U = 1 m/sec, H = 100 m
                dil_fac = np.array(dilution_factor_sectorwise_all_distances, dtype=object).T
                if self.config['Y'] == 0 and self.config['Z'] == 0 and self.config['measurement_height'] == 100 \
                        and self.config['instantaneous_release_bq_list'] == [31536000] and \
                        self.config['release_height'] == 100 and not self.config['like_to_scale_with_mean_speed:']:
                    self.test_single_plume_glc_hukkoo(dil_fac)
                f.write('\n\n')
                f.write(
                    'Instantaneous release: DILUTION FACTORS(s/m^3) per unit release for 6 stability categories at '
                    'various distances:\n')
                f.write(df_dilution_factor.__repr__())
                if self.config['pickle_it']:
                    df_dilution_factor.to_pickle('%s/dilution_factors_%s.p'% (self.config['logdir_name'], self.config['input_file_name']))
                f.write('\n')
                if self.config['Y'] == 0 and self.config['Z'] == 0:
                    f.write(
                        '\nNote (a): Dilution factor is obtained with Y= 0, Z=0 that projects maximum concentration along '
                        'the projection of plume centerline on the ground.')
                else:
                    f.write(
                        '\nNote (a): Dilution factor is obtained with Y= {}, Z={}.'.format(self.config['Y'],
                                                                                           self.config['Z']))

                f.write('\n\n')
                df_max_dilution_factor = pd.DataFrame(
                    np.max(np.array(dilution_factor_sectorwise_all_distances, dtype=object), axis=1).T,
                    columns=['Max Dilution Factor'], index=[str(dist) + 'm' for dist in self.downwind_distances])
                ndx = np.array(dilution_factor_sectorwise_all_distances, dtype=object).argmax(axis=1)
                max_by_stab_cat = [stab_cat_name[i] for i in ndx]
                df_max_dilution_factor['Stability Category'] = max_by_stab_cat
                f.write(df_max_dilution_factor.__repr__())
                f.write('\n\n')
            if self.config['have_dilution_factor']:
                max_dilution_factor = self.config['list_max_dilution_factor']
                df_max_dilution_factor = pd.DataFrame(max_dilution_factor.items(),
                                                      columns=['Distance', 'Max Dilution Factor'])
                df_max_dilution_factor = df_max_dilution_factor.set_index('Distance')
                f.write(df_max_dilution_factor.__repr__())
                f.write('\n\n')

        if self.config['long_term_release']:
            if not self.config['have_dilution_factor']:
                # conservative approach
                if not self.config['have_met_data']:
                    df_dilution_factor = pd.DataFrame(
                        np.array(dilution_factor_sectorwise_all_distances, dtype=object).T,
                        columns=[str(dist) + 'm' for dist in self.downwind_distances],
                        index=stab_cat_name)
                    f.write('\n\n')
                    f.write(
                        'Long-term release (Conservative approach): DILUTION FACTORS(s/m^3) per unit release for 6 '
                        'stability categories at various distances:\n')

                    f.write(df_dilution_factor.__repr__())
                    f.write('\n\n')
                    if self.config['pickle_it']:
                        df_dilution_factor.to_pickle('%s/dilution_factors_%s.p' % (self.config['logdir_name'], self.config['input_file_name']))
                    df_max_dilution_factor = pd.DataFrame(
                        np.max(np.array(dilution_factor_sectorwise_all_distances, dtype=object), axis=1).T,
                        columns=['Max Dilution Factor'], index=[str(dist) + 'm' for dist in self.downwind_distances])
                    ndx = np.array(dilution_factor_sectorwise_all_distances, dtype=object).argmax(axis=1)
                    max_by_stab_cat = [stab_cat_name[i] for i in ndx]
                    df_max_dilution_factor['Stability Category'] = max_by_stab_cat
                    f.write(df_max_dilution_factor.__repr__())
                    f.write('\n')
                    f.write('Due to unavailability of metereological data, as a conservative approach sector average'
                            'plume concentration corresponding for six stability categories are computed for wind '
                            'speed = 1 m/s. Finally, stability category which gives highest concentration is printed. '
                            'Note that: These conditions result in minimal dilution and high plume centerline doses, '
                            ' but also very narrow plumes. (NUREG 1140, US NRC)')
                    f.write('\n\n')

                else:
                    df_dilution_factor = pd.DataFrame(
                        np.array(dilution_factor_sectorwise_all_distances, dtype=object).T,
                        columns=[str(dist) + 'm' for dist in self.downwind_distances],
                        index=sectors)
                    f.write('\n')
                    f.write('Long term release: DILUTION FACTORS(s/m^3) per unit release for 16 wind sectors at '
                            'various distances:\n')
                    f.write(df_dilution_factor.__repr__())
                    f.write('\n\n')
                    # pickle it
                    if self.config['pickle_it']:
                        df_dilution_factor.to_pickle('%s/dilution_factors_%s.p' % (self.config['logdir_name'], self.config['input_file_name']))
                    f.write('\n\n')
                    if self.config['Z'] == 0:
                        f.write(
                            '\nNote (a): Dilution factor is obtained with Z=0 that projects maximum concentration along '
                            'the projection of plume centerline on the ground.')
                    else:
                        f.write(
                            '\nNote (a): Dilution factor is obtained with Z={}.'.format(self.config['Z']))

                    df_max_dilution_factor = pd.DataFrame(
                        np.array(dilution_factor_sectorwise_all_distances, dtype=object).max(axis=1).T,
                        columns=['Max Dilution Factor'], index=[str(dist) + 'm' for dist in self.downwind_distances])
                    ndx = np.array(dilution_factor_sectorwise_all_distances, dtype=object).argmax(axis=1)
                    max_sector = [sectors[i] for i in ndx]
                    df_max_dilution_factor['Sector'] = max_sector
                    f.write(df_max_dilution_factor.__repr__())
                    f.write('\n\n')

            if self.config['have_dilution_factor']:
                max_dilution_factor = self.config['list_max_dilution_factor']
                df_max_dilution_factor = pd.DataFrame(max_dilution_factor.items(),
                                                      columns=['Distance', 'Max Dilution Factor'])
                df_max_dilution_factor = df_max_dilution_factor.set_index('Distance')
                f.write(df_max_dilution_factor.__repr__())
                f.write('\n\n')

        if self.config['run_dose_computation']:
            ### PRINT DCFs ##################
            f.write('\n\n')
            f.write('########## Dose Conversion Factors:##############')
            f.write('\n')
            f.write(
                'Ref 1: External Exposure to Radionuclides in Air, Water and Soil, Federal Guidance Report No. 15. (2019)')
            f.write('\n')
            f.write(
                'Ref 2: International Atomic Energy Agency, Radiation Protection and Safety of Radiation Sources: International Basic Safety Standards. IAEA Safety Standards Series No. GSR Part 3 Vienna (2014)')
            f.write('\n')

            half_lives, lamb_rads = self.find_half_life_and_decay_const_radionuclides()
            df_hl = pd.DataFrame(half_lives, index=self.rads_list, columns=['Half life (in second)'])
            f.write('\n\n')
            f.write('Half life of radionuclides:')
            f.write('\n')
            f.write(df_hl.__repr__())
            print('DCFSSSSS:', DCFs)
            B = np.array(DCFs, dtype=object).T
            for ndx, A in enumerate(B):
                df_dcf = pd.DataFrame(A.T, index=self.age_group,
                                      columns=['Inhalation (Sv/Bq)', 'Ground-shine (Sv-m^2/Bq-second)',
                                               'Submersion (Sv-m^3/Bq-second)', 'Ingestion (Sv/Bq)'])
                f.write('\n\n')
                f.write('Dose Conversion Factors used in the dose computations for radionuclide %s \n' % (
                    str(self.rads_list[ndx])))
                f.write('\n')

                f.write(df_dcf.__repr__())
            f.write('\nNote (a): For ground-shine and submersion dose, progeny corrected DCF values are printed first '
                    'followed by uncorrected DCF values. If there is no progeny for consideration, value remains the same.')
            f.write(
                '\nNote (b): If radionuclide is Tritium (H-3), the values under ingestion indicate dose coefficients for HTO and OBT, respectively.')
            f.write('\n\n')

            ### PRINT DOSE DATA ##################
            if self.config['long_term_release']:
                f.write('################## RESULTS OF DOSE COMPUTATION (in mSv/year) #########################')
            else:
                f.write('################## RESULTS OF DOSE COMPUTATION (in mSv) #########################')
            doses = np.array(DOSES, dtype=object)
            num_dosetype = 3
            dosetype = ['Inhalation', 'Ground-shine', 'Submersion']
            pathtype = ['vegetables', 'meat', 'milk']
            columns = [str(dist) + 'm' for dist in self.downwind_distances]
            rownames = []
            for i in self.rads_list:
                for j in dosetype:
                    for k in self.age_group:
                        index = str(i) + '-' + str(j) + '-' + str('age') + str(k)
                        rownames.append(index)

            df_dose = pd.DataFrame(
                doses[:][:][:].T.reshape(num_dosetype * len(self.rads_list) * len(self.age_group),
                                         len(self.downwind_distances)), columns=columns, index=rownames)

            ### RESULTS OF DOSE COMPUTATION
            f.write('\n\n')
            f.write(
                '########### Inhalation, Ground-shine and Submersion ###########')
            p_igsub = {}
            for dst, _ in enumerate(DOSES):
                for ag, each in enumerate(_):
                    f.write('\n\n')
                    f.write('Downwind Distance: {} metre; Age: {}'.format(self.downwind_distances[dst], self.age_group[ag]))
                    f.write('\n\n')
                    igsub = pd.DataFrame(each, dtype=float, index=['Inhalation', 'Ground-shine', 'Submersion'],
                                         columns=self.rads_list).T
                    igsub['total'] = igsub.loc[:, :].sum(axis=1)
                    f.write(igsub.__repr__())
                    # pickle it
                    if self.config['pickle_it']:
                        igsub.to_pickle('%s/inh_gs_sub_doses_d_%s_age_%s_%s.p' % (self.config['logdir_name'], self.downwind_distances[dst], self.age_group[ag], self.config['input_file_name']))

                    if self.config['downwind_distances'][dst] == self.config['plant_boundary']:
                        p_igsub[ag] = igsub.iloc[:, :3]
            # print('INGESTION_DOSES:', np.array(INGESTION_DOSES))
            ids = np.array(INGESTION_DOSES, dtype=float).reshape(len(self.downwind_distances), len(self.age_group),
                                                                 len(self.rads_list), len(pathtype))

            f.write('\n\n')
            f.write(
                '########### Ingestion Dose ###########')

            # checking for radionuclides whether transfer factor for terrestrial food is available
            lambda_s_per_d_list, lambda_w_per_d_list, f_v1_list, f_v2_list, Fm_Milk_d_per_L_list, Ff_Meat_d_per_kg_list, notransfer_factor_rad = self.fv_list_ecerman_ingestion(
                master_file='library/Dose_ecerman_final.xlsx', sheet_name='eco_param')
            if not pd.DataFrame(notransfer_factor_rad).empty:
                f.write('\n\n')
                f.write("WARNING: ingestion dose values (milk/meat) not computed for the following elements: {}"
                        " as the ELEMENT SPECIFIC TRANSFER FACTORS FOR TERRESTRIAL FOODS is/are not available".format(
                    notransfer_factor_rad))
                f.write('\n')

            total_ingestion = {}
            # print('INGESTION_DOSES:', np.array(INGESTION_DOSES).shape)
            for dst, _ in enumerate(ids):
                for ag, each in enumerate(_):
                    f.write('\n\n')
                    f.write('Downwind Distance: {} metre; Age: {}'.format(self.downwind_distances[dst], self.age_group[ag]))
                    f.write('\n\n')
                    # print("each", each)
                    ings = pd.DataFrame(each.T, dtype=float, index=['Veg', 'Milk', 'Meat'], columns=self.rads_list).T
                    ings = self.zeroing_ingestion(ings, notransfer_factor_rad)
                    ings['total'] = ings.loc[:, :].sum(axis=1)
                    f.write(ings.__repr__())
                    # pickle it
                    if self.config['pickle_it']:
                        ings.to_pickle('%s/ings_d_%s_age_%s_%s.p' % (self.config['logdir_name'], self.downwind_distances[dst], self.age_group[ag], self.config['input_file_name']))

                    if self.config['downwind_distances'][dst] == self.config['plant_boundary']:
                        # ids = np.array(ing).reshape(2, 2, 2, 3)
                        # ping = ings['total']
                        total_ingestion[ag] = ings['total']

            total_inhal_gs_submers_dose_per_rad_list = []
            for age in self.age_group:
                string = 'age' + str(age)
                df = df_dose[df_dose.index.str.contains(r'\b%s\b' % string)]
                rownames = []
                for i in self.rads_list:
                    for j in dosetype:
                        for k in self.age_group:
                            index = str(i) + '-' + str(j)
                        rownames.append(index)
                df.index = rownames

                for rad in self.rads_list:
                    # total inhalation_dose age wise
                    string = str(rad) + '-' + 'Inhalation'
                    df_boundary = df[pboundary]
                    df_rad = df_boundary[df_boundary.index.str.contains(r'\b%s\b' % (string))]
                    # total ground shine age wise
                    string = str(rad) + '-' + 'Ground-shine'
                    total_inhal_gs_submers_dose_per_rad_list.append(np.array(df_rad).flatten()[0])
                    df_rad = df_boundary[df_boundary.index.str.contains(r'\b%s\b' % string)]
                    # total submersion age wise
                    string = str(rad) + '-' + 'Submersion'
                    total_inhal_gs_submers_dose_per_rad_list.append(np.array(df_rad).flatten()[0])
                    df_rad = df_boundary[df_boundary.index.str.contains(r'\b%s\b' % string)]
                    total_inhal_gs_submers_dose_per_rad_list.append(np.array(df_rad).flatten()[0])

                # print plume shine doses for specific radionuclide for 6 STABILITY CATEGORIES at specified distances
                energies, emmission_prob, en_dict, neglected_energies = self.gamma_energy_abundaces(
                    master_file="library/Dose_ecerman_final.xlsx",
                    sheet_name="gamma_energy_radionuclide")

                if (self.config['run_plume_shine_dose'] and age == self.age_group[-1] and self.config['single_plume']) \
                        or (self.config['long_term_release'] and age == self.age_group[-1]
                            and self.config['run_plume_shine_dose'] and not self.config['have_met_data']):

                    f.write('\n\n')
                    if self.config['long_term_release'] and not self.config['have_met_data']:
                        f.write('########### PLUME SHINE DOSE (For Six Stability Categories) (Sector-averaged Plume) '
                                '###########')
                    else:
                        f.write('########### PLUME SHINE DOSE (For Six Stability Categories) (Instanteneous Release) '
                                '###########')
                    f.write('\n')
                    if len(neglected_energies) > 0:
                        f.write(
                            "following gamma energies are neglected for plume shine dose computation (in dict format): {}".format(
                                neglected_energies))
                    for ndx, each in enumerate(PLUME_DOSES):
                        f.write('\n\n')
                        f.write('Downwind Distance: {} metre'.format(self.downwind_distances[ndx]))
                        df_plume_dose = pd.DataFrame(each.reshape(len(self.rads_list), 6), dtype=float,
                                                     index=self.rads_list,
                                                     columns=['A', 'B', 'C', 'D', 'E', 'F'])
                        f.write('\n\n')
                        f.write(df_plume_dose.__repr__())
                        # pickle it
                        if self.config['pickle_it']:
                            df_plume_dose.to_pickle('%s/ings_d_%s_%s.p' % (self.config['logdir_name'],
                            self.downwind_distances[ndx], self.config['input_file_name']))

                # print plume shine doses for specific radionuclide for 16 sectors at specified distances
                elif self.config['run_plume_shine_dose'] and age == self.age_group[-1] \
                        and self.config['long_term_release'] and self.config['have_met_data']:
                    f.write('\n\n')
                    f.write('########### PLUME SHINE DOSE for radionuclides (Long Term Release Scenario) ###########')
                    f.write('\n')
                    if len(neglected_energies) > 0:
                        f.write(
                            "following gamma energies are neglected for plume shine dose computation (in dict format): {}".format(
                                neglected_energies))
                    # averaging plume shine dose over all years of meteorological data
                    for ndx, each in enumerate(PLUME_DOSES):
                        f.write('\n\n')
                        f.write('Downwind Distance: {} metre'.format(self.downwind_distances[ndx]))
                        each = np.array(each).sum(axis=0)
                        df_plume_dose = pd.DataFrame(each.reshape(len(self.rads_list), 16), dtype=float,
                                                     index=self.rads_list,
                                                     columns=sectors)
                        df_plume_dose['Maximum'] = df_plume_dose.loc[:, :].max(axis=1)
                        df_plume_dose['Maximum Sector'] = df_plume_dose.loc[:, :].idxmax(axis=1)
                        f.write('\n\n')
                        f.write(df_plume_dose.__repr__())
                        # pickle it
                        if self.config['pickle_it']:
                            df_plume_dose.to_pickle('%s/ings_d_%s_%s.p' % (self.config['logdir_name'],
                            self.downwind_distances[ndx], self.config['input_file_name']))

                else:
                    pass

            f.write('\n\n')
            f.write('###########RESULTS OF TOTAL DOSE COMPUTATION AT PLANT BOUNDARY###########')

            for _, age in enumerate(self.age_group):
                pigsub = p_igsub[_]
                ping = total_ingestion[_]
                df_all = pd.concat((pigsub, ping), axis=1)
                df_all.columns = df_all.columns.str.replace('total', 'Ingestion (Sum of all pathways)')
                df_all['Total'] = df_all.loc[:, :].sum(axis=1)
                f.write('\n\n')
                if self.config['long_term_release']:
                    f.write('Calculated total dose values (mSv/year) for age %s at specified plant boundary (%s metre)' % (
                        age, self.config['plant_boundary']))
                else:
                    f.write('Calculated total dose values (mSv) for age %s at specified plant boundary (%s metre)' % (
                        age, self.config['plant_boundary']))
                f.write('\n\n')
                f.write(df_all.__repr__())
            f.write('\n\n')
        f.close()
        return

    def agewise_dcfs_inh_gs_submersion(self, age):
        """
            Extract age-wise dose conversion factors (DCFs) for inhalation, ground shine, and submersion.

            Parameters:
                age (int): Age group for which to extract DCFs.

            Returns:
                list: List of DCFs for inhalation, ground shine, and submersion, along with ingestion DCFs.

            Notes:
                This method extracts age-wise dose conversion factors (DCFs) for inhalation, ground shine, and submersion
                from specified master files and sheets.

        """
        DCFs = []
        logging.getLogger("DCF extraction:").info(
            "Age: {age} \n\n".format(age=age))
        # inhalation
        list_dcf_inhalation = self.inhalation_dcf_list(master_file="library/RadioToxicityMaster.xls",
                                                       sheet_name="Inhalation CED Sv per Bq Public",
                                                       age=age)
        # ground shine
        list_dcf_ecerman_ground_shine = self.dcf_list_ecerman_ground_shine_include_progeny(
            master_file="library/Dose_ecerman_final.xlsx", sheet_name="surface_dose", age=age,
            consider_progeny=True)
        # submersion
        list_dcf_ecerman_submersion = self.dcf_list_ecerman_submersion_include_progeny(
            master_file="library/Dose_ecerman_final.xlsx",
            sheet_name="submersion_dose", age=age)
        # ingestion
        list_dcf_ingestion = self.dcf_list_ingestion(master_file="library/Dose_ecerman_final.xlsx",
                                                     sheet_name="ingestion_gsr3", age=age)
        DCFs.append([list_dcf_inhalation, list_dcf_ecerman_ground_shine, list_dcf_ecerman_submersion, list_dcf_ingestion])
        return DCFs

    def agewise_dose_inh_gs_submersion(self, X1, age):
        """
            Calculate age-wise doses for inhalation, ground shine, and submersion at a specified distance.

            Parameters:
                X1 (float): Distance at which doses are calculated.
                age (int): Age group for which doses are calculated.

            Returns:
                list: List of age-wise doses for inhalation, ground shine, and submersion.

            Notes:
                This method calculates age-wise doses for inhalation, ground shine, and submersion at a specified distance.
                It uses the maximum dilution factor for the distance to compute the doses.

        """
        DOSES_AGEWISE = []
        logging.getLogger("Dose calculation:").info(
            "Age: {age}, Distance: {dist}\n\n".format(age=age, dist=X1))
        max_chi_by_q = self.max_dil_fac_all_dist.loc[str(X1)][0]
        # inhalation
        inhalation_dose_values = self.inhalation_dose(X1, age=age, max_dilutfac_for_distance_secperm3=max_chi_by_q)
        inhalation_dose_values = np.array(inhalation_dose_values).flatten()
        # ground shine
        ground_shine_dose_values = self.ground_shine_dose(X1, age=age, max_dilutfac_for_distance_secperm3=max_chi_by_q)
        ground_shine_dose_values = np.array(ground_shine_dose_values).flatten()
        # submersion
        list_submersion_dose = self.submersion_dose(X1, age=age, max_dilutfac_for_distance_secperm3=max_chi_by_q)
        list_submersion_dose = np.array(list_submersion_dose).flatten()

        DOSES_AGEWISE.append([inhalation_dose_values, ground_shine_dose_values, list_submersion_dose])
        return DOSES_AGEWISE

    def agewise_ingestion_dose(self, X1, age):
        """
            Calculate age-wise ingestion dose at a specified distance.

            Parameters:
                X1 (float): Distance at which doses are calculated.
                age (int): Age in years for which ingestion dose is calculated.

            Returns:
                list: List containing the ingestion dose for the specified age.

            Notes:
                This method calculates age-wise ingestion dose at a specified distance based on the age of the individual.
                It uses the maximum dilution factor for the distance to compute the dose.

        """
        max_chi_by_q = self.max_dil_fac_all_dist.loc[str(X1)][0]
        INGESTION_DOSES_AGEWISE = []
        if age > 17:
            ingestion_dose_values = self.ingestion_dose(X1, receiver='adult', max_dilutfac_for_distance_secperm3=max_chi_by_q)
        elif age == 1:
            ingestion_dose_values = self.ingestion_dose(X1, receiver='infant', max_dilutfac_for_distance_secperm3=max_chi_by_q)
        else:
            pass
        INGESTION_DOSES_AGEWISE.append(ingestion_dose_values)
        return INGESTION_DOSES_AGEWISE

    # plume-shine dose
    def agewise_plume_shine_dose(self, X1, age):
        PLUME_DOSES = []
        if self.config['run_plume_shine_dose']:
            # age plume shine dose is age-invariant so will be performed for last entry in age_group
            # last_age = [each for each in self.age_group if each == self.age_group[-1]][0]
            if age == self.age_group[-1]:
                # list_gamma_energy_abundaces = self.gamma_energy_abundaces(
                #     master_file="library/Dose_ecerman_final.xlsx",
                #    sheet_name="gamma_energy_radionuclide")
                print(
                    "age-invariant plume shine dose is being calculated for downwind distance of {}".format(
                        X1))
                plumeshine_dose_values = self.plumeshine_dose(spatial_distance=X1)
                PLUME_DOSES.append(plumeshine_dose_values)
        return PLUME_DOSES

    # plume-shine dose
    def all_dist_agewise_plume_shine_dose(self):
        """
            Calculate age-wise plume-shine dose at a specified distance.

            Parameters:
                X1 (float): Distance at which doses are calculated.
                age (int): Age in years for which plume-shine dose is calculated.

            Returns:
                list: List containing the plume-shine dose for the specified age.

            Notes:
                This method calculates age-wise plume-shine dose at a specified distance based on the age of the individual.
                Plume-shine dose calculation is performed only if the configuration allows it (config['run_plume_shine_dose'] is True).
                Plume-shine dose is age-invariant and is calculated for the last entry in the age_group list.

        """
        if self.config['run_plume_shine_dose']:
            # age plume shine dose is age-invariant so will be performed for last entry in age_group
            # last_age = [each for each in self.age_group if each == self.age_group[-1]][0]
            PLUME_DOSES = []
            for X1 in self.downwind_distances:
                print("age-invariant plume shine dose is being calculated for downwind distance of {}".format(
                            X1))
                plumeshine_dose_values = self.plumeshine_dose(spatial_distance=X1)
                PLUME_DOSES.append(plumeshine_dose_values)
        else:
            PLUME_DOSES = None
        return PLUME_DOSES

    def parallel_allage_dose_inh_gs_submersion(self):
        """
            Calculate doses for all combinations of downwind distances and age groups in parallel.

            Returns:
                list: List containing dose values for all combinations of downwind distances and age groups.

            Notes:
                This method parallelizes the calculation of doses for all combinations of downwind distances and age groups.
                It uses the Parallel function from the joblib library with the maximum number of jobs set to the number of available CPU cores.
                Each combination of downwind distance and age group is processed independently, improving overall computation speed.
        """
        doses = Parallel(n_jobs=-1, verbose=100)(
            delayed(self.agewise_dose_inh_gs_submersion)(dst, ag) for dst, ag in
            list(itertools.product(self.downwind_distances, self.age_group)))
        return list(doses)

    def parallel_allage_dose_ingestion(self):
        """
            Calculate ingestion doses for all combinations of downwind distances and age groups in parallel.

            Returns:
                list: List containing ingestion dose values for all combinations of downwind distances and age groups.

            Notes:
                This method parallelizes the calculation of ingestion doses for all combinations of downwind distances and age groups.
                It uses the Parallel function from the joblib library with the maximum number of jobs set to the number of available CPU cores.
                Each combination of downwind distance and age group is processed independently, improving overall computation speed.
        """
        ing_doses = Parallel(n_jobs=-1, verbose=100)(
            delayed(self.agewise_ingestion_dose)(dst, ag) for dst, ag in
            list(itertools.product(self.downwind_distances, self.age_group)))
        return list(ing_doses)

    def parallel_ps_dose(self):
        """
            Calculate plume-shine doses for all combinations of downwind distances and age groups in parallel.

            Returns:
                list: List containing plume-shine dose values for all combinations of downwind distances and age groups.

            Notes:
                This method parallelizes the calculation of plume-shine doses for all combinations of downwind distances and age groups.
                It uses the Parallel function from the joblib library with the maximum number of jobs set to the number of available CPU cores.
                Each combination of downwind distance and age group is processed independently, improving overall computation speed.
        """
        ps_doses = Parallel(n_jobs=-1, verbose=100)(
            delayed(self.agewise_plume_shine_dose)(dst, ag) for dst, ag in
            list(itertools.product(self.downwind_distances, self.age_group)))
        return list(ps_doses)

    def parallel_allage_dcfs_inh_gs_submersion(self):
        """
            Calculate dose conversion factors (DCFs) for inhalation, ground shine, and submersion for all age groups in parallel.

            Returns:
                list: List containing dose conversion factors for inhalation, ground shine, and submersion for all age groups.

            Notes:
                This method parallelizes the calculation of dose conversion factors (DCFs) for inhalation, ground shine, and submersion
                for all age groups. It uses the Parallel function from the joblib library with the maximum number of jobs set to the
                number of available CPU cores. Each age group is processed independently, improving overall computation speed.
        """
        dcfs = Parallel(n_jobs=-1, verbose=100)(
            delayed(self.agewise_dcfs_inh_gs_submersion)(ag) for ag in self.age_group)
        return list(dcfs)

    def dil_fac_all_sectors_all_dist(self, X1):
        """
            Calculate the dilution factor for all sectors at a given downwind distance.

            Args:
                X1 (float): The downwind distance at which the dilution factor is calculated.

            Returns:
                dict or None: A dictionary containing the dilution factor for each sector at the specified downwind distance.
                    If the configuration specifies that dilution factors are not available, returns None.

            Notes:
                This method calculates the dilution factor for each sector at a given downwind distance. The dilution factor
                represents the reduction in concentration of a substance as it disperses in the atmosphere. If dilution factors
                are not available based on the configuration, it returns None. This method is useful for understanding the
                dispersion pattern of substances released into the environment.
        """
        if not self.config['have_dilution_factor']:
            dilution_factor_sectorwise = self.dilution_per_sector(X1)
            # print('dilution_factor_sectorwise:', dilution_factor_sectorwise)
            # if self.config['long_term_release'] and self.config['have_met_data']:
                # if self.config['plot_dilution_factor']:
                # plot bar plot of dilution factor comprising 16 sectors
            #    self.plot_dilution_factor(X1)
            #    self.max_dilution_factor = max(dilution_factor_sectorwise)
        else:
            dilution_factor_sectorwise = None
        return dilution_factor_sectorwise

    def parallel_dilfac_alldists(self):
        """
            Calculate the dilution factor for all sectors at multiple downwind distances in parallel.

            Returns:
                list: A list containing the dilution factor for each downwind distance. Each element in the list is a dictionary
                    containing the dilution factor for each sector at the corresponding downwind distance.

            Notes:
                This method calculates the dilution factor for each sector at multiple downwind distances in parallel using
                parallel processing. It iterates over each downwind distance and calculates the dilution factors for all sectors
                at that distance. The result is a list where each element corresponds to a downwind distance and contains a
                dictionary with the dilution factor for each sector at that distance.
        """
        dilfac = Parallel(n_jobs=-1, verbose=100)(
            delayed(self.dil_fac_all_sectors_all_dist)(d) for d in self.downwind_distances)
        return list(dilfac)

    def dose_calculation_script(self):
        """
        Run script for dose calculation.

        This method executes the entire process of dose calculation, including preprocessing, missing data correction,
        calm correction factor calculation, parallel computation of dilution factors, dose computation for inhalation,
        ground shine, submersion, and ingestion, as well as plume shine dose computation if enabled.

        Returns:
            tuple or dict: Depending on the configuration, the method returns either a tuple or a dictionary containing
                the calculated results. If the configuration specifies not to run dose computation, it returns only the
                dilution factor for all sectors at all distances. If dose computation is enabled, it returns a tuple
                containing various dose-related data, including DCFs, dilution factors, doses for inhalation, ground shine,
                submersion, and ingestion. If plume shine dose computation is also enabled, it includes plume shine doses
                in the return value.

        Notes:
            - This method first checks the configuration settings to determine which parts of the dose calculation process
              should be executed.
            - It runs different preprocessing steps based on the configuration settings.
            - If parallel computation is enabled, it uses parallel processing to calculate dilution factors, inhalation,
              ground shine, submersion, and ingestion doses simultaneously.
            - It saves the computed data using pickle if configured to do so.
            - Finally, it returns the computed results as specified by the configuration.
        """
        if self.config['single_plume']:
            self.TJFD_ALL_MISSING_CORR = self.synthetic_TJFD_for_single_plume()

        if self.config['long_term_release']:
            self.file_preprocessing()
            self.met_data_to_tjfd()
            self.missing_correction()
            self.calm_correction_factor_calc()

        if self.config['long_term_release'] and self.config['have_met_data'] and not self.config['have_dilution_factor']:
            start_operation_time = self.config['start_operation_time']
            end_operation_time = self.config['end_operation_time']
            self.operation_hours_per_day = np.absolute(start_operation_time - end_operation_time)

        header = ''
        for ndx, each in enumerate(self.downwind_distances):
            if ndx == 0:
                header += 'dist(m):' + '%1.0f' % each
            else:
                header += '%11.0f' % each
        #################### PARALLEL IMPLEMENTATION ##################################################################
        # SHAPES OF DOSES INGESTION_DOSES, DCFS, PLUME_DOSES, dilution_factor_sectorwise_all_distances
        # PLUME_DOSES = len(self.downwind_distances), len(self.rads_list), len(self.sheets), 1, 1, 16
        # DCFs.shape:(len(self.age_group), 4, 1)
        # DOSES.shape: (len(self.downwind_distances), len(self.age_group), 3, 1)
        # INGESTION_DOSES.shape: (len(self.downwind_distances), len(self.age_group), 1, 3)
        # dilution_factor_sectorwise_all_distances.shape(len(self.downwind_distances), 16) for
        # long_term_no_met otherwise len(self.downwind_distances),6
        # DILUTION_FACTOR_ALL_SECTORS_ALL_DISTANCES
        dilution_factor_sectorwise_all_distances = self.parallel_dilfac_alldists()
        self.max_dil_fac_all_dist = self.get_max_dilution_factor(dilution_factor_sectorwise_all_distances)
        if self.config['run_dose_computation']:
            DCFs = self.parallel_allage_dcfs_inh_gs_submersion()
            DCFs = np.array(DCFs, dtype=object).reshape(len(self.age_group), 4, len(self.rads_list))
            DOSES = self.parallel_allage_dose_inh_gs_submersion()
            DOSES = np.array(DOSES, dtype=object).reshape(len(self.downwind_distances), len(self.age_group), 3, len(self.rads_list))
            INGESTION_DOSES = self.parallel_allage_dose_ingestion()
            INGESTION_DOSES = np.array(INGESTION_DOSES, dtype=object).reshape(len(self.downwind_distances), len(self.age_group), len(self.rads_list), 3)
            PLUME_DOSES = self.all_dist_agewise_plume_shine_dose()

        # PICKLE_IT
        if self.config['pickle_it']:
            #  dilution_factor_sectorwise_all_distances.shape(len(dist), len(sectors)/len(stab_cat))
            with open("%s/Dil_Fac_%s.p" % (self.config['logdir_name'], self.config['input_file_name']), "wb") as output_file:
                cPickle.dump(dilution_factor_sectorwise_all_distances, output_file)
            if self.config['run_dose_computation']:
                # DOSES.shape: (len(self.downwind_distances), len(self.age_group), 3, 1)
                with open("%s/I_G_S_Doses_%s.p" % (self.config['logdir_name'], self.config['input_file_name']), "wb") as output_file:
                    cPickle.dump(DOSES, output_file)
                # INGESTION_DOSES.shape: (len(self.downwind_distances), len(self.age_group), 1, 3)
                with open("%s/ING_Doses_%s.p" % (self.config['logdir_name'], self.config['input_file_name']), "wb") as output_file:
                    cPickle.dump(INGESTION_DOSES, output_file)
                # PLUME_DOSES = len(self.downwind_distances), len(self.rads_list), len(self.sheets), 1, 1, 16
                if self.config['run_plume_shine_dose']:
                    with open("%s/PS_Doses_%s.p" % (self.config['logdir_name'], self.config['input_file_name']), "wb") as output_file:
                        cPickle.dump(PLUME_DOSES, output_file)
        ############################## NOT PARALLEL ####################################################################
        # dilution_factor_sectorwise_all_distances, DCFs, DOSES, PLUME_DOSES, INGESTION_DOSES = self.simple_loop_output()
        if not self.config['run_dose_computation']:
            return dilution_factor_sectorwise_all_distances
        if self.config['run_dose_computation']:
            if self.config['run_plume_shine_dose']:
                return DCFs, dilution_factor_sectorwise_all_distances, DOSES, INGESTION_DOSES, PLUME_DOSES
            else:
                return DCFs, dilution_factor_sectorwise_all_distances, DOSES, INGESTION_DOSES

    def simple_loop_output(self):
        """
            Perform dose calculation using a simple loop approach.

            This method computes dilution factors, DCFs, doses for inhalation, ground shine, submersion, and ingestion using
            a simple loop approach without parallelization.

            Returns:
                tuple: A tuple containing the computed results, including dilution factors for all sectors at all distances,
                    DCFs for different exposure pathways and age groups, doses for inhalation, ground shine, submersion, and
                    ingestion for different downwind distances, plume shine doses if enabled, and ingestion doses for
                    different age groups and distances.
        """
        dilution_factor_sectorwise_all_distances = []
        DOSES = []
        PLUME_DOSES = []
        INGESTION_DOSES = []
        for X1 in np.array(self.downwind_distances):
            if not self.config['have_dilution_factor']:
                # distance based inhalation
                dilution_factor_sectorwise = self.dilution_per_sector(X1)
                # print('dilution_factor_sectorwise:', dilution_factor_sectorwise)
                if self.config['long_term_release'] and self.config['have_met_data']:
                    # if self.config['plot_dilution_factor']:
                    # plot bar plot of dilution factor comprising 16 sectors
                    self.plot_dilution_factor(X1)
                    self.max_dilution_factor = max(dilution_factor_sectorwise)
                dilution_factor_sectorwise_all_distances.append(dilution_factor_sectorwise)

            if self.config['run_dose_computation']:
                # DCFs depend on age
                DCFs = []
                DOSES_AGEWISE = []
                INGESTION_DOSES_AGEWISE = []
                for age in self.age_group:
                    logging.getLogger("Dose calculation:").info(
                        "Age: {age}, Distance: {dist}\n\n".format(age=age, dist=X1))
                    # inhalation
                    list_dcf_inhalation = self.inhalation_dcf_list(master_file="library/RadioToxicityMaster.xls",
                                                                   sheet_name="Inhalation CED Sv per Bq Public",
                                                                   age=age)
                    inhalation_dose_values = self.inhalation_dose(X1, age=age)
                    inhalation_dose_values = np.array(inhalation_dose_values).flatten()
                    # ground shine
                    list_dcf_ecerman_ground_shine = self.dcf_list_ecerman_ground_shine_include_progeny(
                        master_file="library/Dose_ecerman_final.xlsx", sheet_name="surface_dose", age=age,
                        consider_progeny=True)
                    ground_shine_dose_values = self.ground_shine_dose(X1, age=age)
                    ground_shine_dose_values = np.array(ground_shine_dose_values).flatten()
                    # submersion
                    list_dcf_ecerman_submersion = self.dcf_list_ecerman_submersion_include_progeny(
                        master_file="library/Dose_ecerman_final.xlsx",
                        sheet_name="submersion_dose", age=age)
                    # print('list_dcf_ecerman_submersion with progenies contribution:', list_dcf_ecerman_submersion)
                    list_submersion_dose = self.submersion_dose(X1, age=age)
                    list_submersion_dose = np.array(list_submersion_dose).flatten()
                    # ingestion
                    list_dcf_ingestion = self.dcf_list_ingestion(master_file="library/Dose_ecerman_final.xlsx",
                                                                 sheet_name="ingestion_gsr3", age=age)

                    if age > 17:
                        ingestion_dose_values = self.ingestion_dose(X1, receiver='adult')
                    elif age == 1:
                        ingestion_dose_values = self.ingestion_dose(X1, receiver='infant')
                    else:
                        pass

                    # plume-shine dose
                    if self.config['run_plume_shine_dose']:
                        # age plume shine dose is age-invariant so will be performed for last entry in age_group
                        # last_age = [each for each in self.age_group if each == self.age_group[-1]][0]
                        if age == self.age_group[-1]:
                            list_gamma_energy_abundaces = self.gamma_energy_abundaces(
                                master_file="library/Dose_ecerman_final.xlsx",
                                sheet_name="gamma_energy_radionuclide")
                            print(
                                "age-invariant plume shine dose is being calculated for downwind distance of {}".format(
                                    X1))
                            plumeshine_dose_values = self.plumeshine_dose(spatial_distance=X1)
                            PLUME_DOSES.append(plumeshine_dose_values)

                    DCFs.append(
                        [list_dcf_inhalation, list_dcf_ecerman_ground_shine, list_dcf_ecerman_submersion,
                         list_dcf_ingestion])

                    DOSES_AGEWISE.append([inhalation_dose_values, ground_shine_dose_values, list_submersion_dose])

                    INGESTION_DOSES_AGEWISE.append(ingestion_dose_values)
                # self.max_dilution_factor = None
                DOSES.append(DOSES_AGEWISE)
                INGESTION_DOSES.append(INGESTION_DOSES_AGEWISE)
        return dilution_factor_sectorwise_all_distances, DCFs, DOSES, PLUME_DOSES, INGESTION_DOSES
