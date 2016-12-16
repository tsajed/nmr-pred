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
    ///NmrExperiment exp = new NmrExperiment();
    callMatlab(folder);
  }

  public static void runPrediction(File folder) {
    try {
    //LinearRegression model = (LinearRegression) weka.core.SerializationHelper.read("models/regression.model_3d");
      //RandomForest modesl = (RandomForest) weka.core.SerializationHelper.read("models/classification.model_1");

      ArrayList<NmrStructure> structures = NmrExperiment.createNmrStructures(folder);
      ArrayList<String> hmdb_ids = NmrExperiment.getStructures(structures, folder);

      Instances isTrainingSet = (Instances) weka.core.SerializationHelper.read("models/train_classification_6");
      RandomForest model = new RandomForest();
      model.buildClassifier(isTrainingSet);

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

  public static void callMatlab(Instances test, ArrayList<String> hmdb_ids) throws 
      MatlabConnectionException, MatlabInvocationException {
        
    MatlabProxyFactoryOptions options =
            new MatlabProxyFactoryOptions.Builder()
                .setUsePreviouslyControlledSession(true)
                .build();
    MatlabProxyFactory factory = new MatlabProxyFactory(options);
    MatlabProxy proxy = factory.getProxy();

    proxy.eval("cd ..");
    proxy.eval("cd matlab");
    // Create variables like sys and inter
    proxy.feval("create_nmr1H_plot", "");
    //proxy.eval("exit()");

    proxy.disconnect();
  }

}