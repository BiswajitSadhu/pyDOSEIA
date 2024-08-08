import sys
import os
import random
import numpy as np
import argparse
# import torch
import logging
import yaml
import subprocess
# import torch
import logging
from outputfunc import *
# import odfpy
import ast

sys.path.append(".")

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inpdir", type=str, help="input file path")
    parser.add_argument("--library", type=str, help="path for library")
    parser.add_argument("--device", choices=["cuda", "cpu"], default="cpu")
    parser.add_argument("--logdir", required=True, type=str, help="info.log will be generated")
    # parser.add_argument("--seed", type=int, help="random seed for weight initialization")
    # parser.add_argument("--output_file_name", required=True, type=str, help="output file name", default="pydoseia.txt")
    directory = parser.parse_args().logdir
    if not os.path.exists(directory):
        os.makedirs(directory)
    return parser.parse_args()


def parse_config(configfile: str):
    with open(configfile, "r") as f:
        data = yaml.load(f, Loader=yaml.FullLoader)
        # yaml.load(f, Loader=yaml.FullLoader)
        return convert_none_to_str(data)


def convert_none_to_str(data):
    if isinstance(data, list):
        data[:] = [convert_none_to_str(i) for i in data]
    elif isinstance(data, dict):
        for k, v in data.items():
            data[k] = convert_none_to_str(v)
    return 'None' if data is None else data

def test1_single_plume_glc_hukkoo(config, dilution_factor_sectorwise_all_distances):
    dil_fac = np.array(dilution_factor_sectorwise_all_distances, dtype=object).sum(axis=2).T
    dil_fac = np.array(dil_fac, dtype=float).T
    hukkoo_glc_100m_single_plume = np.array([[2.2840147833303766e-14, 5.590981589039797e-07, 1.5015465916440028e-05,
                                              1.8534496171567465e-05, 9.310282356034462e-06, 3.6973162422584492e-06,
                                              9.083901472613404e-07, 4.64985088543347e-07, 5.809360379016258e-08],
                                             [5.841151371898279e-22, 1.881938226103609e-09, 8.545950204498517e-07,
                                              1.2296547461666177e-05, 1.689222515044389e-05, 1.3596757373581253e-05,
                                              6.948090247268007e-06, 4.718541203801888e-06, 1.2717701082737338e-06],
                                             [7.956707894464337e-42, 1.091680781377853e-14, 2.634808537371747e-09,
                                              1.504493947251421e-06, 7.156250389638748e-06, 1.2763487559346726e-05,
                                              1.1755542375765323e-05, 9.483132536764242e-06, 3.5384718119021803e-06],
                                             [1.9691262342948032e-107, 1.7554619140561087e-32, 2.322138749483794e-18,
                                              1.645624503785968e-10, 3.9826750086653584e-08, 8.73472028784478e-07,
                                              4.620161828899931e-06, 6.349041184912812e-06, 6.810849638824183e-06],
                                             [6.22620696059002e-181, 9.192930646978123e-57, 1.5133024136726587e-31,
                                              1.0210504603823029e-16, 6.745515294647785e-12, 5.643704797281375e-09,
                                              4.787367464209608e-07, 1.3403801979441503e-06, 4.725940039623769e-06],
                                             [0.0, 2.5748671022239412e-139, 1.0642330561676144e-74,
                                              9.189248563624478e-36, 5.821816484706907e-23, 3.877805752304435e-15,
                                              5.70587584195314e-10, 8.896730306387878e-09, 5.399181471824446e-07]],
                                            dtype=float).T

    dict_hukkoo_test_glc = {}
    test_dist = [100, 200, 300, 500, 700, 1000, 1600, 2000, 4000]
    for d, glc in zip(test_dist, hukkoo_glc_100m_single_plume):
        dict_hukkoo_test_glc[d] = glc
    dict_dil_fac = {}
    for d, glc in zip([100, 200, 300, 500, 700, 1000, 1600, 2000, 4000], dil_fac):
        dict_dil_fac[d] = glc
    # run_test = False
    for k in dict_dil_fac.keys():
        if k in test_dist:
            # run_test = True
            assert np.allclose(dict_dil_fac[k], dict_hukkoo_test_glc[k])
            print("TEST 1 passed (Verified with Hukkoo and Bapat results)!!!"
                  " on computing ground level concentration for height 100 m "
                  "with release rate 1 Bq/S, mean speed 1 m/S for downwind distance {} m".format(k))
            logging.getLogger("TEST1").info(
                "config file name: test1_single_conc_hukkoo.yaml; content: {config_desc}".format(config_desc=config))
            logging.getLogger("TEST1").info(
                "TEST 1 passed (Verified with Hukkoo and Bapat results)!!!"
                " on computing ground level concentration for height 100 m "
                "with release rate 1 Bq/S, mean speed 1 m/S for downwind distance {} m".format(k))
    return


def test2_inhal_gs_submersion(config, doses):
    print("TEST 2 DETAILS: DILUTION FACTOR AVAILABLE; INPUT FILE NAME: test_inhal_gs.yaml. This checks that the code "
          "produce correct inhalation, ground-shine and submersion dose age 18; expected inhalation, ground-shine "
          "and submersion dose ")
    print('config for TEST 2:', config)
    logging.getLogger("TEST2").info(
        "config file name test_inhal_gs.yaml; content: {config_desc}".format(config_desc=config))
    logging.getLogger("TEST2_details").info(
        "TEST 2 DETAILS: DILUTION FACTOR AVAILABLE; INPUT FILE NAME: test_inhal_gs.yaml. This checks that the code "
        "produce correct inhalation, ground-shine and submersion dose age 18; expected inhalation, ground-shine "
        "and submersion dose")

    igs_1 = [np.array([7.14817841e-09, 1.19136307e-08, 2.72311558e-08, 7.14817841e-08,
                       1.22540201e-08, 1.87214196e-08]),
             np.array([9.99870151e-06, 4.33399982e-06, 5.84129433e-06, 5.25050336e-08,
                       1.43664677e-08, 3.18601397e-06]),
             np.array([5.75063659e-10, 3.54622590e-10, 2.87148454e-10, 1.89004256e-12,
                       8.85598035e-11, 1.29214688e-10])]

    # exp_igsub_18 = pd.DataFrame(igs_18, dtype=float, index=['Inhalation', 'Ground-shine', 'Submersion'],
    #                           columns=['Co-60', 'Cs-134', 'Eu-154', 'Sr-90', 'I-131', 'Cs-137']).T

    # age 18
    igs_18 = [np.array([1.02116834e-08, 2.04233669e-08, 5.41219222e-08, 1.63386935e-07,
                        7.55664575e-09, 3.98255654e-08]),
              np.array([7.89641042e-06, 3.37916548e-06, 4.55736627e-06, 4.30066356e-08,
                        1.08526877e-08, 2.48452454e-06]),
              np.array([4.52383412e-10, 2.69129792e-10, 2.20824445e-10, 1.54500436e-12,
                        6.47905056e-11, 9.77575020e-11])]
    # exp_igsub_1 = pd.DataFrame(igs_1, dtype=float, index=['Inhalation', 'Ground-shine', 'Submersion'],
    #                            columns=['Co-60', 'Cs-134', 'Eu-154', 'Sr-90', 'I-131', 'Cs-137']).T

    for dst, _ in enumerate(doses):
        for ag, each in enumerate(_):
            # igsub = pd.DataFrame(each, dtype=float, index=['Inhalation', 'Ground-shine', 'Submersion'],
            #                     columns=['Co-60', 'Cs-134', 'Eu-154', 'Sr-90', 'I-131', 'Cs-137']).T
            if ag == 0:
                assert np.allclose(igs_1, each)
                # print("Test passed (Infant dose Computation)! The inhalation, ground-shine and submersion modules "
                #      "are working as expected for infant age group.")
                logging.getLogger("TEST2_results").info(
                    "Test passed (Infant dose Computation)! The inhalation, ground-shine and submersion modules "
                    "are working as expected for infant age group.")

            if ag == 1:
                assert np.allclose(igs_18, each)
                # print("Test passed (Adult dose computation)! The inhalation, ground-shine and submersion modules "
                #      "are working as expected for adult age group.")
                logging.getLogger("TEST2_results").info(
                    "Test passed (Adult dose Computation)! The inhalation, ground-shine and submersion modules "
                    "are working as expected for adult age group.")

    return



def test3_ingestion_iaea_iv11(config, ingestion_doses):
    print("TEST 3 DETAILS: INGESTION DOSE COMPARISION; INPUT FILE NAME: iv_11.yaml. This checks that the code "
          "produce correct ingestion dose AT 1 KM for age 1 and 18 for I-131 RELEASE at 60 m ELEVATION HEIGHT."
          "and submersion dose ")
    print('config for TEST 3:', config)
    logging.getLogger("TEST3").info(
        "config file name test_ingestion.yaml; content: {config_desc}".format(config_desc=config))
    logging.getLogger("TEST3_details").info(
        "TEST 3 DETAILS: INGESTION DOSE COMPARISION; INPUT FILE NAME: iv_11.yaml. This checks that the code "
        "produce correct ingestion dose AT 1 KM for age 1 and 18 for I-131 RELEASE at 60 m ELEVATION HEIGHT."
        "and submersion dose ")

    ings_c = np.array([[[[2.313066e-05, 0.00015630349896783523, 1.5194067250317405e-05]],

              [[7.705844677039068e-06, 1.591941343058461e-05, 4.641784309634899e-06]]]])

    assert np.allclose(ingestion_doses, ings_c)
    print("Test passed (Infant dose Computation)! The ingestion modules are working as" 
        " expected for infant age group.")
    logging.getLogger("TEST3_results").info(
                    "Test passed (Ingestion dose Computation)! The ingestion modules "
                    "are working as expected for infant age group.")

    return

def test_main():
    """
            Main function to parse arguments, configure logging, parse configuration file, perform dose calculation, and generate output.
    """
    ########## TEST 1: Test Dilution factor segment ####################################################################
    args = parse_arguments()
    # seed_torch(seed=int(args.seed))
    log_file_name = args.config_file.split('.yaml')[0]
    logging.basicConfig(filename='%s_info.log' % log_file_name, level=logging.INFO, filemode="w+")
    logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))
    # logging.getLogger("main").info(f"cuda available: {torch.cuda.is_available()}")
    config = parse_config('test/test1_single_conc_hukkoo.yaml')
    # metfunc = MetFunc(config=config, device=args.device, logdir=args.logdir)
    output_func = OutputFunc(config=config, device=args.device, logdir=args.logdir)
    # all_results = output_func.dose_calculation_script()

    if not config['run_dose_computation']:
        dilution_factor_sectorwise_all_distances = output_func.dose_calculation_script()
        output_func.output_to_txt(dilution_factor_sectorwise_all_distances, filename=args.output_file_name,
                                  DCFs=None, DOSES=None, INGESTION_DOSES=None, PLUME_DOSES=None)

    if config['run_dose_computation']:
        if config['run_plume_shine_dose']:
            DCFs, dilution_factor_sectorwise_all_distances, DOSES, INGESTION_DOSES, PLUME_DOSES = output_func.dose_calculation_script()
            test1_single_plume_glc_hukkoo(config, dilution_factor_sectorwise_all_distances)

        else:
            DCFs, dilution_factor_sectorwise_all_distances, DOSES, INGESTION_DOSES = output_func.dose_calculation_script()
            test1_single_plume_glc_hukkoo(config, dilution_factor_sectorwise_all_distances)

    logging.getLogger("main").info("path of input file: {config_desc}".format(config_desc=args.config_file))
    logging.getLogger("main").info(
        "output file name: {output_file_name}".format(output_file_name=args.output_file_name))

    ########## TEST 2: Inhalation, Ground Shine Dose ###################################################################
    args = parse_arguments()
    # seed_torch(seed=int(args.seed))
    log_file_name = args.config_file.split('.yaml')[0]
    logging.basicConfig(filename='%s_info.log' % log_file_name, level=logging.INFO, filemode="w+")
    logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))
    # logging.getLogger("main").info(f"cuda available: {torch.cuda.is_available()}")
    config = parse_config('test/test_inhal_gs.yaml')
    # metfunc = MetFunc(config=config, device=args.device, logdir=args.logdir)
    output_func = OutputFunc(config=config, device=args.device, logdir=args.logdir)
    # all_results = output_func.dose_calculation_script()

    if not config['run_dose_computation']:
        dilution_factor_sectorwise_all_distances = output_func.dose_calculation_script()
        output_func.output_to_txt(dilution_factor_sectorwise_all_distances, filename=args.output_file_name,
                                  DCFs=None, DOSES=None, INGESTION_DOSES=None, PLUME_DOSES=None)

    if config['run_dose_computation']:
        if config['run_plume_shine_dose']:
            DCFs, dilution_factor_sectorwise_all_distances, DOSES, INGESTION_DOSES, PLUME_DOSES = output_func.dose_calculation_script()
            test2_inhal_gs_submersion(config, DOSES)

        else:
            DCFs, dilution_factor_sectorwise_all_distances, DOSES, INGESTION_DOSES = output_func.dose_calculation_script()
            test2_inhal_gs_submersion(config, DOSES)

    ####### TEST 3: Ingestion DOSE with IAEA GSR Excercise IV-11 #######################################################
    args = parse_arguments()
    # seed_torch(seed=int(args.seed))
    log_file_name = args.config_file.split('.yaml')[0]
    logging.basicConfig(filename='%s_info.log' % log_file_name, level=logging.INFO, filemode="w+")
    logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))
    # logging.getLogger("main").info(f"cuda available: {torch.cuda.is_available()}")
    config = parse_config('test/iv_11.yaml')
    # metfunc = MetFunc(config=config, device=args.device, logdir=args.logdir)
    output_func = OutputFunc(config=config, device=args.device, logdir=args.logdir)
    # all_results = output_func.dose_calculation_script()

    if not config['run_dose_computation']:
        dilution_factor_sectorwise_all_distances = output_func.dose_calculation_script()
        output_func.output_to_txt(dilution_factor_sectorwise_all_distances, filename=args.output_file_name,
                                  DCFs=None, DOSES=None, INGESTION_DOSES=None, PLUME_DOSES=None)

    if config['run_dose_computation']:
        if config['run_plume_shine_dose']:
            DCFs, dilution_factor_sectorwise_all_distances, DOSES, INGESTION_DOSES, PLUME_DOSES = output_func.dose_calculation_script()
            test3_ingestion_iaea_iv11(config, INGESTION_DOSES)

        else:
            DCFs, dilution_factor_sectorwise_all_distances, DOSES, INGESTION_DOSES = output_func.dose_calculation_script()
            test3_ingestion_iaea_iv11(config, INGESTION_DOSES)

if __name__ == "__main__":
    test_main()
