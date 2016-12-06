import java.lang.String;
import java.io.*;
import java.util.*;
import java.util.regex.*;
import java.nio.file.*;
import java.nio.charset.StandardCharsets;
import javax.swing.JFrame;
import java.awt.*;
import java.text.DecimalFormat;

import weka.core.*;
import weka.classifiers.functions.*;
import weka.classifiers.functions.supportVector.*;
import weka.classifiers.Evaluation;
import weka.classifiers.evaluation.Prediction;
import weka.classifiers.trees.*;
import weka.classifiers.*;
import weka.attributeSelection.*;
import weka.classifiers.meta.*;

import org.math.plot.*;
import org.math.io.*;


public class NmrPred {

  public static void main(String[] argv) {
  	File folder = new File("dataset/");
    try {
      //LinearRegression model = (LinearRegression) weka.core.SerializationHelper.read("models/regression.model_3d");
      Instances isTrainingSet = (Instances) weka.core.SerializationHelper.read("models/train_classification_2");
      //runRegression(isTrainingSet, true); 
      runClassifier(isTrainingSet, true);
    }
    catch (Exception e) {
      e.printStackTrace();
      //Instances trainingSet = buildTrainingClassification(folder);
      Instances trainingSet = buildTrainingRegression(folder);
      runRegression(trainingSet, false);
      //runClassifier(trainingSet, false);
    }    
  }

  static void runClassifier(Instances isTrainingSet, boolean read) {
    try {
      RandomForest d_tree_model = new RandomForest();
      //Bagging d_tree_model = new Bagging();

      //Setting Parameters
      // d_tree_model.setLearningRate(0.1);
      // d_tree_model.setMomentum(0.2);
      // d_tree_model.setTrainingTime(500);
      //d_tree_model.setHiddenLayers("100");
      // SMO d_tree_model = new SMO();
      //d_tree_model.buildClassifier(isTrainingSet);
      if (!read) {
        weka.core.SerializationHelper.write("models/classification.model_1", d_tree_model);
        weka.core.SerializationHelper.write("models/train_classification_2", isTrainingSet);
      }


      for (int i = 25; i < 116; i++) {
        int j = 25;
        //for (int j = 0; j < 20; j++) {
          //System.out.println(isTrainingSet.instance(j).value(i));
        isTrainingSet.deleteAttributeAt(j);
        //}
      }

      for (int i = 0; i < 21; i++) {
        int j = 0;
        //for (int j = 0; j < 20; j++) {
          //System.out.println(isTrainingSet.instance(j).value(i));
        isTrainingSet.deleteAttributeAt(j);
        //}
      }

      for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 200; j++) {
          System.out.println(isTrainingSet.instance(j).value(i));
        }
      }

      //isTrainingSet = performFeatureExtraction(isTrainingSet);
      d_tree_model.buildClassifier(isTrainingSet);
      Evaluation eTest = new Evaluation(isTrainingSet);
     // eTest.evaluateModel(d_tree_model, isTrainingSet);
      Random rand = new Random(1);
      eTest.crossValidateModel(d_tree_model, isTrainingSet, 10, rand);
      String strSummary = eTest.toSummaryString();
      System.out.println(strSummary);
      ArrayList<Prediction> predictions = eTest.predictions();

      double true_values[] = new double[predictions.size()];
      double predicted_values[] = new double[predictions.size()];

      double error = 0;
      for (int i = 0; i < predictions.size(); i++) {
        true_values[i] = predictions.get(i).actual();
        predicted_values[i] = predictions.get(i).predicted();
        if (Math.abs(predictions.get(i).predicted() - predictions.get(i).actual()) < 8) {
          //error = Math.abs(true_values[i] - predicted_values[i]) + error;
        }
        error = Math.abs(true_values[i] - predicted_values[i]) + error;
      }

      error = error / predictions.size();
      System.out.println(error);

      Plot2DPanel plot = new Plot2DPanel();

      plot.addScatterPlot("Linear Scatter Plot", Color.RED, true_values, predicted_values);
      plot.addLinePlot("True Regression Plot", Color.BLUE, true_values, true_values);

      // put the PlotPanel in a JFrame, as a JPanel
      JFrame frame = new JFrame("Panel");
      frame.setContentPane(plot);
      frame.setVisible(true);

    }
    catch (Exception e) {
      e.printStackTrace();
    }
  }

  static void runRegression(Instances isTrainingSet, boolean read)  {
    try {
      //model.buildClassifier(isTrainingSet);
      //LinearRegression model = new LinearRegression();
      MultilayerPerceptron model = new MultilayerPerceptron();
      //Setting Parameters
      model.setLearningRate(0.1);
      model.setMomentum(0.2);
      model.setTrainingTime(2000);
      model.setHiddenLayers("116,116,116,116");
      if (!read) {
        weka.core.SerializationHelper.write("models/regression.model_1", model);
        weka.core.SerializationHelper.write("models/train_regression_1", isTrainingSet);
      }

      // for (int i = 29; i < 116; i++) {
      //   int j = 29;
      //   isTrainingSet.deleteAttributeAt(j);
      // }


      //isTrainingSet = performFeatureExtraction(isTrainingSet);
      model.buildClassifier(isTrainingSet);
      Evaluation eTest = new Evaluation(isTrainingSet);
      //eTest.evaluateModel(model, isTrainingSet);
      Random rand = new Random(1);
      eTest.crossValidateModel(model, isTrainingSet, 5, rand);
      String strSummary = eTest.toSummaryString();
      System.out.println(strSummary);
      ArrayList<Prediction> predictions = eTest.predictions();

      double true_values[] = new double[predictions.size()];
      double predicted_values[] = new double[predictions.size()];

      for (int i = 0; i < predictions.size(); i++) {
        if (Math.abs(predictions.get(i).predicted() - predictions.get(i).actual()) > 5) {

        }
        true_values[i] = predictions.get(i).actual();
        predicted_values[i] = predictions.get(i).predicted();
      }

      Plot2DPanel plot = new Plot2DPanel();

      plot.addScatterPlot("Linear Scatter Plot", Color.RED, true_values, predicted_values);
      plot.addLinePlot("True Regression Plot", Color.BLUE, true_values, true_values);

      // put the PlotPanel in a JFrame, as a JPanel
      JFrame frame = new JFrame("Panel");
      frame.setContentPane(plot);
      frame.setVisible(true);

    }
    catch (Exception e) {
      e.printStackTrace();
    }
  }

  static Instances performFeatureExtraction(Instances data) {
    //PrincipalComponents pcaEvaluator = new PrincipalComponents();
    //WrapperSubsetEval evaluator = new WrapperSubsetEval();
    System.out.println(data.numAttributes());
    CfsSubsetEval evaluator = new CfsSubsetEval();
    int k = data.numAttributes();
    J48 classifier = new J48();
    //evaluator.setClassifier(classifier);

    // Sets the amount of variance to account for when retaining principal
    // components.
    //pcaEvaluator.setVarianceCovered(1);
    // Sets maximum number of attributes to include in transformed attribute
    // names.
    //pcaEvaluator.setMaximumAttributeNames(-1);

    // Scaled X such that the variance of each feature is 1.
    //boolean scale = Utils.getFlag('s', args);
    //pcaEvaluator.setCenterData(true);

      //pcaEvaluator.setCenterData(false);
    // Ranking the attributes.
    BestFirst ranker = new BestFirst();
    // Specify the number of attributes to select from the ranked list.
    //ranker.setNumToSelect(k - 1);

    try {
      AttributeSelection selector = new AttributeSelection();
      selector.setSearch(ranker);
      selector.setEvaluator(evaluator);
      selector.SelectAttributes(data);

      // Transform data into eigenvector basis.
      Instances transformedData = selector.reduceDimensionality(data);
      System.out.println(transformedData.numAttributes());
      return transformedData;
    } catch (Exception e) {
      e.printStackTrace();
      return data;
    }

  }

  static Instances buildTrainingRegression(File folder) {
    int feature_factor = 4;

    ArrayList<NmrStructure> nmr_structures = getChemicalShifts(folder);
    getStructures(nmr_structures, folder);
    for (NmrStructure nmr_str : nmr_structures) {
      try {
        System.out.println(nmr_str.hmdb_id);
        nmr_str.atomic_descriptors = GetCDKDescriptors.getAtomicDescriptor(nmr_str.structure_sdf, "");
        nmr_str.findNearestAtomToHydrogens(GetCDKDescriptors.getNearestAtoms(nmr_str.structure_sdf));
      }
      catch (Exception e)
      {
        e.printStackTrace();
      }
    }
    ArrayList<Double[]> values = nmr_structures.get(0).atomic_descriptors;
    ArrayList<Attribute> attributes = new ArrayList<Attribute>();
    for (int i = 0; i < (feature_factor * values.size())+1; i++) {
      attributes.add(new Attribute(String.valueOf(i)));
    }
    FastVector wekaAttributes = new FastVector(feature_factor*values.size()+1);
    for(Attribute a : attributes) {
      wekaAttributes.addElement(a);
    }
    Instances isTrainingSet = new Instances("Rel", wekaAttributes, 500);
    isTrainingSet.setClassIndex(feature_factor*values.size());
    System.out.println(feature_factor*values.size());

    /* i = carbon positions from nmr shift text file
       j = feature index 
       atomic_descriptor.get(j)[i] gets jth feature (descriptor) and ith position in hydrogen_positions
       namely ith Carbon atom as numbered by molfile and text file

       k is how many feature factors we are taking as descriptors
       If we take 2 nearest atoms, so all 29 features of those 2 nearest atoms are included. k = 3
    */

    for (NmrStructure nmr_str : nmr_structures) {
      for (int i = 0; i < nmr_str.hydrogen_positions.size(); i++) {
        Instance iExample = new DenseInstance(feature_factor*values.size() + 1);
        for (int j = 0; j < nmr_str.atomic_descriptors.size(); j++) {
          iExample.setValue((Attribute)wekaAttributes.elementAt(j), 
                            nmr_str.atomic_descriptors.get(j)[Integer.valueOf(nmr_str.hydrogen_positions.get(i))]);
        }
        for (int k = 0; k < feature_factor - 1; k++) {
          for (int j = 0; j < nmr_str.atomic_descriptors.size(); j++) {
            iExample.setValue((Attribute)wekaAttributes.elementAt(j + (k+1)*values.size()), 
                              nmr_str.atomic_descriptors.get(j)[Integer.valueOf(nmr_str.nearest_atoms.get(i).get(k))]);
          }
        }
        iExample.setValue((Attribute)wekaAttributes.elementAt(feature_factor*values.size()), nmr_str.chemical_shifts.get(i));

        isTrainingSet.add(iExample);
      }
    }
    return isTrainingSet;
  }

  static Instances buildTrainingClassification(File folder) {
    int feature_factor = 4;
    ArrayList<NmrStructure> nmr_structures = getChemicalShifts(folder);
    getStructures(nmr_structures, folder);
    for (NmrStructure nmr_str : nmr_structures) {
      try {
        System.out.println(nmr_str.hmdb_id);
        nmr_str.atomic_descriptors = GetCDKDescriptors.getAtomicDescriptor(nmr_str.structure_sdf, "");
        nmr_str.findNearestAtomToHydrogens(GetCDKDescriptors.getNearestAtoms(nmr_str.structure_sdf));
      }
      catch (Exception e)
      {
        e.printStackTrace();
      }
    }
    //nmr_structures.get(0).printDescriptors();
    ArrayList<Double[]> values = nmr_structures.get(0).atomic_descriptors;
    ArrayList<Attribute> attributes = new ArrayList<Attribute>();
    for (int i = 0; i < feature_factor* values.size(); i++) {
      attributes.add(new Attribute(String.valueOf(i)));
    }

    FastVector fv = new FastVector(100);
    DecimalFormat df = new DecimalFormat("#.#");
    for (float i = 1; i <= 100; i++) {
        fv.addElement(df.format(i/10));
        System.out.println(df.format(i/10));
    } 
    attributes.add(new Attribute("Class", fv));
    FastVector wekaAttributes = new FastVector(feature_factor * values.size()+1);

    for(Attribute a : attributes) {
      wekaAttributes.addElement(a);
    }
    Instances isTrainingSet = new Instances("Rel", wekaAttributes, 500);
    isTrainingSet.setClassIndex(feature_factor*values.size());

    for (NmrStructure nmr_str : nmr_structures) {
      for (int i = 0; i < nmr_str.hydrogen_positions.size(); i++) {
        Instance iExample = new DenseInstance(feature_factor*values.size() + 1);
        for (int j = 0; j < nmr_str.atomic_descriptors.size(); j++) {
          iExample.setValue((Attribute)wekaAttributes.elementAt(j), 
                            nmr_str.atomic_descriptors.get(j)[Integer.valueOf(nmr_str.hydrogen_positions.get(i))]);
        }
        for (int k = 0; k < feature_factor - 1; k++) {
          for (int j = 0; j < nmr_str.atomic_descriptors.size(); j++) {
            iExample.setValue((Attribute)wekaAttributes.elementAt(j + (k+1)*values.size()), 
                              nmr_str.atomic_descriptors.get(j)[Integer.valueOf(nmr_str.nearest_atoms.get(i).get(k))]);
          }
        }

        iExample.setValue((Attribute)wekaAttributes.elementAt(feature_factor*values.size()), nmr_str.c_shift_classes.get(i));
        isTrainingSet.add(iExample);
      }
    }
    return isTrainingSet;
  }

  static ArrayList<NmrStructure> getChemicalShifts(File folder) {
    ArrayList<NmrStructure> structures = new ArrayList<NmrStructure>();

    for (File fileEntry : folder.listFiles()) {
      ArrayList<String> carbon_position = new ArrayList<String>();
      ArrayList<Float> chemical_shift = new ArrayList<Float>();

      String name = fileEntry.getName();
      String pattern = ".txt$";
      String[] file_names = name.split("\\.");
      // Create a Pattern object
      Pattern r = Pattern.compile(pattern);

      // Now create matcher object.
      Matcher m = r.matcher(name);
      if (m.find( )) {
        readChemicalShifts(fileEntry, carbon_position, chemical_shift);
        NmrStructure structure = new NmrStructure(carbon_position, chemical_shift, "", file_names[0]);
        structures.add(structure);
      }
    }

    return structures;
  }

  static void readChemicalShifts(File file, ArrayList<String> c_pos, ArrayList<Float> c_shift) {
    BufferedReader reader = null;
    try {
      reader = new BufferedReader(new FileReader(file));
      String text = null;
      boolean parse_shifts = false;
      while ((text = reader.readLine()) != null) {
        String[] items = text.split("\t");
        if (parse_shifts) {
          if (text.equals("")) {
            continue;
          }
          if (text.contains("Table")) {
            parse_shifts = false;
            continue;
          }
          System.out.println(file.getName());
          if (items[0].trim().equals("No.") || items[2].contains("M")) {
            continue;
          }

          c_pos.add(items[1].trim());
          if (items[2] != null && !items[2].isEmpty()) {
            c_shift.add(Float.parseFloat(items[2].trim()));
          }
          else {
            c_shift.add(Float.parseFloat(items[3].trim()));
          }
          System.out.println(file.getName());
        }
        if (items.length == 0 || items.length == 1) {
          continue;
        }
        if (items[1].trim().equals("Atom")) {
          parse_shifts = true;
        }
      }
    } catch (FileNotFoundException e) {
        e.printStackTrace();
    } catch (IOException e) {
        e.printStackTrace();
    } finally {
      try {
        if (reader != null) {
          reader.close();
        }
      } catch (IOException e) {
      }
    }
  }

  static void getStructures(ArrayList<NmrStructure> structures, File folder) {
    for (NmrStructure nmr_s: structures) {
      String file = folder.getName() + "/" + nmr_s.hmdb_id + ".sdf";

      try {
        String text = new String(Files.readAllBytes(Paths.get(file)), StandardCharsets.UTF_8);  
        nmr_s.structure_sdf = text;   
      } catch (FileNotFoundException e) {
          e.printStackTrace();
      } catch (IOException e) {
          e.printStackTrace();
      }
    }
  }

}

