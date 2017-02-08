import java.lang.String;
import java.io.*;
import java.util.*;
import java.util.regex.*;
import java.nio.file.*;

import weka.core.*;
import weka.classifiers.functions.*;
import weka.classifiers.functions.supportVector.*;
import weka.classifiers.Evaluation;
import weka.classifiers.evaluation.Prediction;
import weka.classifiers.trees.*;
import weka.classifiers.*;

import matlabcontrol.*;

public class NmrPred {
  public static void main(String[] argv) throws MatlabConnectionException, MatlabInvocationException {
    File folder = new File("test/");
    //runPrediction(folder);
    NmrExperiment exp = new NmrExperiment();
  }

  /* RunPrediction function builds a classifier from a saved Training Set. It then builds a testing instance
   * from a test folder. It then runs predict function to calculate ppm values for all Hydrogen atoms in 
   * for all the test molecules
   */
  public static void runPrediction(File folder) throws MatlabConnectionException, MatlabInvocationException {
    try {
    //LinearRegression model = (LinearRegression) weka.core.SerializationHelper.read("models/regression.model_3d");
      //RandomForest modesl = (RandomForest) weka.core.SerializationHelper.read("models/classification.model_1");

      ArrayList<NmrStructure> structures = NmrExperiment.createNmrStructures(folder);
      ArrayList<String> hmdb_ids = NmrExperiment.getStructures(structures, folder);

      Instances isTrainingSet = (Instances) weka.core.SerializationHelper.read("models/train_classification_6");
      RandomForest model = new RandomForest();
      model.buildClassifier(isTrainingSet);

      Instances test = NmrExperiment.buildTestClassification(folder, structures);

      int i = 0;
      for (NmrStructure nmr_str : structures) {
        for (String h_pos : nmr_str.hydrogen_positions) {
          // Class is a double from 1 to 100 classes. So divide by 10 to get shift
          Double clsLabel = model.classifyInstance(test.instance(i))/10;
          String label = String.valueOf(clsLabel);

          nmr_str.chemical_shifts.add(Float.valueOf(label));
          System.out.println(hmdb_ids.get(i) + "\t" + label);
          i++;
        }
      }

      callMatlab(structures);
    }
    catch (Exception e) {
      e.printStackTrace();
    }
            
  }

  public static void callMatlab(ArrayList<NmrStructure> structures) throws 
      MatlabConnectionException, MatlabInvocationException {

    MatlabProxyFactoryOptions options =
            new MatlabProxyFactoryOptions.Builder()
                .setUsePreviouslyControlledSession(true)
                .build();
    MatlabProxyFactory factory = new MatlabProxyFactory(options);
    MatlabProxy proxy = factory.getProxy();

    proxy.eval("cd ..");
    proxy.eval("cd matlab");

    for (NmrStructure nmr_str : structures) {
      String isotopes = "sys.isotopes = {";
      String ppms = "inter.zeeman.scalar = {";

      // Empty any existing coupling constants from previous runs
      proxy.eval("inter.coupling.scalar = []");
      
      for (int i = 0; i < nmr_str.hydrogen_positions.size(); i++) {
        if (i == nmr_str.hydrogen_positions.size() - 1) {
          isotopes = isotopes + "'1H'}";
          ppms = ppms + String.valueOf(nmr_str.chemical_shifts.get(i)) + "}";
          String coupling = "inter.coupling.scalar{" + String.valueOf(i+1) + "," +
                     String.valueOf(i+1) + "} = 0" ;
          proxy.eval(coupling);
        }
        else {
          isotopes = isotopes + "'1H', ";
          ppms = ppms + String.valueOf(nmr_str.chemical_shifts.get(i)) + ", ";
        }

        for (int j = i + 1; j < nmr_str.hydrogen_positions.size(); j++) {
          String coupling = "inter.coupling.scalar{" + String.valueOf(i+1) + "," + 
                             String.valueOf(j+1) + "} = 1" ;
          proxy.eval(coupling);
        }
      }
      System.out.println(ppms);
      proxy.eval(isotopes);
      proxy.eval(ppms);

      proxy.feval("create_nmr1H_plot");
    }

    proxy.disconnect();
  }

}