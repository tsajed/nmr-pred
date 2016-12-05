import java.util.*;
import java.text.DecimalFormat;

public class NmrStructure {
  ArrayList<String> hydrogen_positions;
  ArrayList<Float> chemical_shifts;
  ArrayList<String> c_shift_classes;
  ArrayList<ArrayList<String>> nearest_atoms;

  String structure_sdf;
  String hmdb_id;
  ArrayList<Double[]> atomic_descriptors;

  public NmrStructure(ArrayList<String> h_pos, ArrayList<Float> c_shifts, String sdf, String id) {
    hydrogen_positions = h_pos;
    chemical_shifts = c_shifts;
    structure_sdf = sdf;
    hmdb_id = id;
    c_shift_classes = new ArrayList<String>();

    assignShiftClasses(chemical_shifts);
  }

  public void assignShiftClasses(ArrayList<Float> c_shifts) {
    DecimalFormat df = new DecimalFormat("#.#");
    for (Float c_shift : c_shifts) {
      c_shift_classes.add(df.format(c_shift));
    }
  }

  // [ [ 0, 4, 2, 1], [ 5, 3, 2, 10] ] => indices of nearest_atoms
  // Only contains hydrogen atoms but they can contain indices of any atoms, sorted by distance
  public void findNearestAtomToHydrogens(ArrayList<ArrayList<String>> n_atom) {
    nearest_atoms = new ArrayList<ArrayList<String>>();
    for (String h_position : hydrogen_positions) {
      ArrayList<String> nearest_pos = n_atom.get(Integer.valueOf(h_position));
      nearest_atoms.add(nearest_pos);
    }
  }

  public void printDescriptors() {
    for (ArrayList<String> h_pos : nearest_atoms) {
      for (String index : h_pos) {
        for (Double[] values : atomic_descriptors) {
          System.out.println(String.valueOf(values[Integer.valueOf(index)]));
        }
        System.out.println("For atom number ---- " + index);
      }
    }
  }
}