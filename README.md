![pydoseia_logo(1)](https://github.com/user-attachments/assets/dfa5fd0f-4d20-4481-8ca8-a90dc9dbe627) 
## pyDOSEIA
**A python-based package for computation of Dose from atmospheric release**

Description: We introduce pyDOSEIA, a robust Python package, designed for radiation risk assessment and dose calculation in scenarios such as nuclear events, radiological accidents, and environmental contamination. Built on advanced computational models, it offers tools for estimating doses from various exposure pathways like inhalation, ingestion, external exposure, and plume shine. Featuring parallel processing and up-to-date dose conversion factors, pyDOSEIA ensures accurate calculations for both short-term and long-term exposures. With a user-friendly interface, it empowers researchers and policymakers in radiation risk assessment and emergency preparedness.

![pydoseia_FLOW_figure1(1)](https://github.com/user-attachments/assets/40e8e66b-af27-4c35-a374-a275d3abf47b)






run Interactive Input Generator using following command:



python auto_input_generator.py

It will prompt questions with hints on how to answer it. Follow it to generate the input file. Once all asked questions are answered by the user, the code automatically initiate computation

If input is already created, user are requested to run auto_input_generator.py and provide the directory name and name of the existing input file name. pyDOSEIA automatically detect the input file and run computation on it.

# optional option: pickle_it == False (default)
With pickle_it:True, pickle files are written using output of dilution factor and dose computations. This is highly useful for analizing the data for research, analysis and machine/deep learning application. 
