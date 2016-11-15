import java.lang.String;
import java.io.*;
import java.util.*;
import java.util.regex.*;
import java.nio.file.*;
import java.nio.charset.StandardCharsets;

import weka.core.*;
import weka.classifiers.functions.*;
import weka.classifiers.Evaluation;


public class NmrPred {

  public static void main(String[] argv) {
  	File folder = new File("dataset/");
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
    LinearRegression model = new LinearRegression();
    try {
      model.buildClassifier(isTrainingSet);
      //Instance test = isTrainingSet.lastInstance();
      //double ppm = model.classifyInstance(test);
      //System.out.println("Predicted ppm = "+ ppm);
      Evaluation eTest = new Evaluation(isTrainingSet);
      // eTest.evaluateModel(model, isTrainingSet);
      // String strSummary = eTest.toSummaryString();
      // System.out.println(strSummary);
      Random rand = new Random(1);
      eTest.crossValidateModel(model, isTrainingSet, 2, rand);
      String strSummary = eTest.toSummaryString();
      System.out.println(strSummary);
      //System.out.println(eTest.toMatrixString());
 
      // Get the confusion matrix (not possible with linear regressor)
      // double[][] cmMatrix = eTest.confusionMatrix();
    }
    catch (Exception e) {
      e.printStackTrace();
    }
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

