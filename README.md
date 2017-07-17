# nmr-pred

Prediction of NMR Spectra from Structure

Predict chemical shifts from structures in smiles or inchi. Use open source software to convert the peaks into a nmr spectograph of peaks to intensity

Use machine learning models to predict NMR Spectra from Structure


# Code Layout and Flow

NmrPred.java 
  - Main Function in Program. You will have to call the main function no matter what you do.
  - Code to Run NmrExperiment to do 10 fold cross validation on a training set. 
  - You can also RunPrediction Tests based on structures in \test folder. 
  - RunPrediction will CallMatlab function to simulate the NMR Spectra based on predicted chemical shifts and structure
