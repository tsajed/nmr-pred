import java.util.*;

public class NmrStructure {
  ArrayList<String> carbon_positions;
  ArrayList<Float> chemical_shifts;
  String structure_sdf;
  String hmdb_id;

  public NmrStructure(ArrayList<String> c_pos, ArrayList<Float> c_shifts, String sdf, String id) {
    carbon_positions = c_pos;
    chemical_shifts = c_shifts;
    structure_sdf = sdf;
    hmdb_id = id;
  }
}