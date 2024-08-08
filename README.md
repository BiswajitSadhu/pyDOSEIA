![pydoseia_logo(1)](https://github.com/user-attachments/assets/dfa5fd0f-4d20-4481-8ca8-a90dc9dbe627) 
## pyDOSEIA
**A python-based package for computation of Dose from atmospheric release**

Description: We introduce pyDOSEIA, a robust Python package, designed for radiation risk assessment and dose calculation in scenarios such as nuclear events, radiological accidents, and environmental contamination. Built on advanced computational models, it offers tools for estimating doses from various exposure pathways like inhalation, ingestion, external exposure, and plume shine. Featuring parallel processing and up-to-date dose conversion factors, pyDOSEIA ensures accurate calculations for both short-term and long-term exposures. With a user-friendly interface, it empowers researchers and policymakers in radiation risk assessment and emergency preparedness.

**A Simple pyDOSEIA flowchart**
![pydoseia_FLOW_figure1(1)](https://github.com/user-attachments/assets/40e8e66b-af27-4c35-a374-a275d3abf47b)

**INGEN (Interactive Input Generator) flowchart**

![pydoseia_test_cases](https://github.com/user-attachments/assets/6be3e18f-e297-43e0-9e0c-7cff780a9da9)

**Installation**

The installation guides for these environments are provided below:

conda create -n pydose python=3.10
conda activate pydose
pip install numpy scipy matplotlib pandas xlrd

**Usage**
run Interactive Input Generator using following command:
python auto_input_generator.py

It will prompt questions with hints on how to answer it. Follow it to generate the input file. Once all asked questions are answered by the user, the code automatically initiate computation

If input is already created, user are requested to run auto_input_generator.py and provide the directory name and name of the existing input file name. pyDOSEIA automatically detect the input file and run computation on it.

optional option: pickle_it == False (default)
With pickle_it:True, pickle files are written using output of dilution factor and dose computations. This is highly useful for analizing the data for research, analysis and machine/deep learning application. 

**Contributing** 
We welcome contributions to the pyDOSEIA project! If you have suggestions, want to improve the code, or have ideas for new features, please create a Pull Request or raise an issue. If you encounter any bugs or have questions about using the package, don't hesitate to contact us. Your feedback and contributions are invaluable in making pyDOSEIA a better tool for the community.

**Contact information**
Copyright(C) 2023 Author and Developer: Dr. Biswajit Sadhu (biswajit.chem001@gmail.com)
