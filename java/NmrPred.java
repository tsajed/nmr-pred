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

public class NmrPred {
  public static void main(String[] argv) {
    File folder = new File("test/");
    //runPrediction(folder);
    NmrExperiment exp = new NmrExperiment();
  }

  public static void runPrediction(File folder) {
    try {
    //LinearRegression model = (LinearRegression) weka.core.SerializationHelper.read("models/regression.model_3d");
      RandomForest model = (RandomForest) weka.core.SerializationHelper.read("models/classification.model_1");

      ArrayList<NmrStructure> structures = NmrExperiment.createNmrStructures(folder);
      ArrayList<String> hmdb_ids = NmrExperiment.getStructures(structures, folder);

      Instances test = NmrExperiment.buildTestClassification(folder, structures);
      for (int i = 0; i < test.numInstances(); i++) {
        Double clsLabel = model.classifyInstance(test.instance(i));
        test.instance(i).setClassValue(clsLabel);
        System.out.println(hmdb_ids.get(i) + "\t" + String.valueOf(clsLabel));
      }
    }
    catch (Exception e) {
      e.printStackTrace();
    }
            
  }
}