# pyDOSEIA_v1: A python-based package for computation of Dose from atmospheric release
# The computation is parallelized using joblib to faster computation and synthetic dataset generation.

run Interactive Input Generator using following command:

python auto_input_generator.py

It will prompt questions with hints on how to answer it. Follow it to generate the input file. Once all asked questions are answered by the user, the code automatically initiate computation

If input is already created, user are requested to run auto_input_generator.py and provide the directory name and name of the existing input file name. pyDOSEIA automatically detect the input file and run computation on it.

# optional option: pickle_it == False (default)
With pickle_it:True, pickle files are written using output of dilution factor and dose computations. This is highly useful for analizing the data for research, analysis and machine/deep learning application. 
