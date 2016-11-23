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
import weka.classifiers.Evaluation;
import weka.classifiers.evaluation.Prediction;
import weka.classifiers.trees.J48;
import weka.classifiers.functions.*;

import org.math.plot.*;
import org.math.io.*;


public class NmrPred {

  public static void main(String[] argv) {
  	File folder = new File("dataset/");
    try {
      //LinearRegression model = (LinearRegression) weka.core.SerializationHelper.read("models/regression.model");
      Instances isTrainingSet = (Instances) weka.core.SerializationHelper.read("models/train_classification");
      //runLinearRegression(model, isTrainingSet, true); 
      runClassifier(isTrainingSet, true);
    }
    catch (Exception e) {
      e.printStackTrace();
      Instances trainingSet = buildTrainingClassification(folder);
      // LinearRegression model = new LinearRegression();
      // runLinearRegression(model, trainingSet, false);
      runClassifier(trainingSet, false);
    }    
  }

  static void runClassifier(Instances isTrainingSet, boolean read) {
    try {
      //J48 d_tree_model = new J48();
      SMO d_tree_model = new SMO();
      //d_tree_model.buildClassifier(isTrainingSet);
      if (!read) {
        weka.core.SerializationHelper.write("models/classification.model", d_tree_model);
        weka.core.SerializationHelper.write("models/train_classification", isTrainingSet);
      }
      d_tree_model.buildClassifier(isTrainingSet);

      Evaluation eTest = new Evaluation(isTrainingSet);
      Random rand = new Random(1);
      eTest.crossValidateModel(d_tree_model, isTrainingSet, 5, rand);
      String strSummary = eTest.toSummaryString();
      System.out.println(strSummary);
      ArrayList<Prediction> predictions = eTest.predictions();

      double true_values[] = new double[predictions.size()];
      double predicted_values[] = new double[predictions.size()];

      double error = 0;
      for (int i = 0; i < predictions.size(); i++) {
        true_values[i] = predictions.get(i).actual();
        predicted_values[i] = predictions.get(i).predicted();
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

  static void runLinearRegression(LinearRegression model, Instances isTrainingSet, boolean read)  {
    try {
      model.buildClassifier(isTrainingSet);
      if (!read) {
        weka.core.SerializationHelper.write("models/regression.model", model);
        weka.core.SerializationHelper.write("models/train_regression", isTrainingSet);
      }

      Evaluation eTest = new Evaluation(isTrainingSet);
      Random rand = new Random(1);
      eTest.crossValidateModel(model, isTrainingSet, 5, rand);
      String strSummary = eTest.toSummaryString();
      System.out.println(strSummary);
      ArrayList<Prediction> predictions = eTest.predictions();

      double true_values[] = new double[predictions.size()];
      double predicted_values[] = new double[predictions.size()];

      for (int i = 0; i < predictions.size(); i++) {
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

  static Instances buildTrainingRegression(File folder) {
    ArrayList<NmrStructure> nmr_structures = getChemicalShifts(folder);
    getStructures(nmr_structures, folder);
    for (NmrStructure nmr_str : nmr_structures) {
      try {
        System.out.println(nmr_str.hmdb_id);
        nmr_str.atomic_descriptors = GetCDKDescriptors.getAtomicDescriptor(nmr_str.structure_sdf, "");
      }
      catch (Exception e)
      {
        e.printStackTrace();
      }
    }
    ArrayList<Double[]> values = nmr_structures.get(0).atomic_descriptors;
    ArrayList<Attribute> attributes = new ArrayList<Attribute>();
    for (int i = 0; i < values.size()+1; i++) {
      attributes.add(new Attribute(String.valueOf(i)));
    }
    FastVector wekaAttributes = new FastVector(values.size()+1);
    for(Attribute a : attributes) {
      wekaAttributes.addElement(a);
    }
    Instances isTrainingSet = new Instances("Rel", wekaAttributes, 500);
    isTrainingSet.setClassIndex(values.size());

    /* i = carbon positions from nmr shift text file
       j = feature index 
       atomic_descriptor.get(j)[i] gets jth feature (descriptor) and ith position in molecule
       namely ith Carbon atom as numbered by molfile and text file
    */

    for (NmrStructure nmr_str : nmr_structures) {
      for (int i = 0; i < nmr_str.carbon_positions.size(); i++) {
        Instance iExample = new DenseInstance(values.size() + 1);
        for (int j = 0; j < nmr_str.atomic_descriptors.size(); j++) {
          System.out.println(String.valueOf(nmr_str.chemical_shifts.get(i)));
          iExample.setValue((Attribute)wekaAttributes.elementAt(j), 
                            nmr_str.atomic_descriptors.get(j)[Integer.valueOf(nmr_str.carbon_positions.get(i))]);
        }
        iExample.setValue((Attribute)wekaAttributes.elementAt(values.size()), nmr_str.chemical_shifts.get(i));

        isTrainingSet.add(iExample);
      }
    }
    return isTrainingSet;
  }

  static Instances buildTrainingClassification(File folder) {
    ArrayList<NmrStructure> nmr_structures = getChemicalShifts(folder);
    getStructures(nmr_structures, folder);
    for (NmrStructure nmr_str : nmr_structures) {
      try {
        System.out.println(nmr_str.hmdb_id);
        nmr_str.atomic_descriptors = GetCDKDescriptors.getAtomicDescriptor(nmr_str.structure_sdf, "");
      }
      catch (Exception e)
      {
        e.printStackTrace();
      }
    }
    ArrayList<Double[]> values = nmr_structures.get(0).atomic_descriptors;
    ArrayList<Attribute> attributes = new ArrayList<Attribute>();
    for (int i = 0; i < values.size(); i++) {
      attributes.add(new Attribute(String.valueOf(i)));
    }

    FastVector fv = new FastVector(100);
    DecimalFormat df = new DecimalFormat("#.#");
    for (float i = 1; i <= 100; i++) {
        fv.addElement(df.format(i/10));
        System.out.println(df.format(i/10));
    } 
    attributes.add(new Attribute("Class", fv));
    FastVector wekaAttributes = new FastVector(values.size()+1);

    for(Attribute a : attributes) {
      wekaAttributes.addElement(a);
    }
    Instances isTrainingSet = new Instances("Rel", wekaAttributes, 500);
    isTrainingSet.setClassIndex(values.size());

    for (NmrStructure nmr_str : nmr_structures) {
      for (int i = 0; i < nmr_str.carbon_positions.size(); i++) {
        Instance iExample = new DenseInstance(values.size() + 1);
        for (int j = 0; j < nmr_str.atomic_descriptors.size(); j++) {
          iExample.setValue((Attribute)wekaAttributes.elementAt(j), 
                            nmr_str.atomic_descriptors.get(j)[Integer.valueOf(nmr_str.carbon_positions.get(i))]);
        }
        System.out.println(String.valueOf(nmr_str.c_shift_classes.get(i)));
        iExample.setValue((Attribute)wekaAttributes.elementAt(values.size()), nmr_str.c_shift_classes.get(i));
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

