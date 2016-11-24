import java.util.*;
import java.text.DecimalFormat;

public class NmrStructure {
  ArrayList<String> carbon_positions;
  ArrayList<Float> chemical_shifts;
  ArrayList<String> c_shift_classes;
  ArrayList<String> nearest_atoms;

  String structure_sdf;
  String hmdb_id;
  ArrayList<Double[]> atomic_descriptors;

  public NmrStructure(ArrayList<String> c_pos, ArrayList<Float> c_shifts, String sdf, String id) {
    carbon_positions = c_pos;
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

  public void findNearestAtomToHydrogens(ArrayList<String> n_atom) {
    nearest_atoms = new ArrayList<String>();
    for (String h_position : carbon_positions) {
      String nearest_pos = n_atom.get(Integer.valueOf(h_position));
      nearest_atoms.add(nearest_pos);
    }
  }
}