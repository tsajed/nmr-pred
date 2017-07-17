# nmr-pred

Prediction of NMR Spectra from Structure

Predict chemical shifts from structures in smiles or inchi. Use open source software to convert the peaks into a nmr spectograph of peaks to intensity

Use machine learning models to predict NMR Spectra from Structure

CDK Library and Descriptors - http://cdk.github.io/cdk/1.4/docs/api/org/openscience/cdk/qsar/descriptors/atomic/package-summary.html

Spinach Library - http://spindynamics.org/Spinach.php

Spinach Documentation - http://spindynamics.org/wiki/index.php?title=Main_Page

WEKA Documentation - http://www.cs.waikato.ac.nz/ml/weka/


# Code Layout and Flow

NmrPred.java 
  - Main Function in Program. You will have to call the main function no matter what you do.
  - Code to Run NmrExperiment to do 10 fold cross validation on a training set. 
  - You can also RunPrediction Tests based on 3D SDF structures in \test folder
  - RunPrediction will CallMatlab function to simulate the NMR Spectra based on predicted chemical shifts and structure

NmrExperiment.java
  - RunClassifier function starts a 10 fold Cross validation based on a training set ( saved or created ). RandomForest algorithm works th e best currently
  - BuildTrainingClassification function builds a training set for classification problem. For every Structure in the dataset folder, it calculates 28 atomic descriptors for every atom in the structure. It also reads chemical shifts for every Hydrogen in the structure. Now, for every Hydrogen it takes the 28 descriptors for that hydrogen and the nearest three hydrogen atoms for it to build a 128 descriptor featureset for classification. There are 100 classes of chemical shifts from 0.1 to 10.0 for every 0.1.
  - BuildTrainingRegression function that builds a training set for regression analysis
  - BuildTestClassification function builds a test set for classification used by RunPrediction in NmrPred.java
  - GetChemicalShifts uses ReadChemicalShift function to read a list of chemical shift text files downloaded from HMDB
  
GetCDKDescriptors.java
  - Uses CDK Package built in Java to calculate atomic properties of atoms, molecules, distances of hydrogen atoms from given 3D SDF structure
  - GetAtomicDescriptor is the main function used to calculate 28 Atomic Descriptors for every atom in a particular structure/molecule. Returns an ArrayList of doubles (descriptors). The Arraylist is the number of atoms in molecule and each molecule will have 28 doubles/descriptors
  - GetNearestAtoms is used to calculate the nearest atoms to any atom in a molecule. This is used to find the three nearest atoms to any Hydrogen atom in a molecule
  - ComputeDescriptorsAtomic calculates atomic descriptors for an Arraylist of Atoms
  - There are molecular descriptors too but they are not currently used for prediction
  
NmrStructure.java
  - A Java class to define an NMR Structure including its chemical shifts of hydrogen atoms, its descriptors for every atom, hydrogen positions, nearest atoms to every atom, hose codes
  - AssignShiftClasses rounds a chemical shift to 0.1 for classification
 
Layout
  - DataSet folder contains the training set of 3D SDF structures and chemical shifts in text files
  - Test Folder contains the test set of 3D SDF structures used for prediction and simulation, not for training or cross-validation
  - Models folder contains saved training classification or regression sets or trained models like RandomForest
  - Matlab Folder contains the only code create_NMR_1H_plot.m in matlab to simulate the NMR spectra based on predicted chemical shifts and known J-Coupling constants. It uses Fourier transformation, apodization etc to create an NMR image.
  - Java folder contains all the prediction algorith, dataset and code in Java
  - Spinach folder contains the Spinach library used for NMR Simulation
  - Docs folder contains relevant papers including my paper for the individual study course
