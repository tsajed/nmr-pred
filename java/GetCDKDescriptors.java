import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Arrays;
import java.util.Vector;
import java.util.Map;

import org.openscience.cdk.ChemFile;
import org.openscience.cdk.aromaticity.CDKHueckelAromaticityDetector;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemFile;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.io.ISimpleChemObjectReader;
import org.openscience.cdk.io.ReaderFactory;
import org.openscience.cdk.qsar.DescriptorEngine;
import org.openscience.cdk.qsar.IDescriptor;
import org.openscience.cdk.qsar.IMolecularDescriptor;
import org.openscience.cdk.qsar.IAtomicDescriptor;
import org.openscience.cdk.qsar.result.DoubleArrayResult;
import org.openscience.cdk.qsar.result.DoubleArrayResultType;
import org.openscience.cdk.qsar.result.DoubleResult;
import org.openscience.cdk.qsar.result.IDescriptorResult;
import org.openscience.cdk.qsar.result.IntegerArrayResult;
import org.openscience.cdk.qsar.result.IntegerArrayResultType;
import org.openscience.cdk.qsar.result.IntegerResult;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;
import org.openscience.cdk.qsar.descriptors.molecular.IPMolecularLearningDescriptor;
import org.openscience.cdk.qsar.descriptors.atomic.IPAtomicLearningDescriptor;
import org.openscience.cdk.smiles.SmilesGenerator;

import javax.vecmath.Point3d;


/**
 * ApplyCDKDescriptors.java
 * Purpose: Calculate CDK descriptors in CSV format from input SDF files
 * For ease of use, the design is completely static, i.e. no member functions
 * Calling the constructor executes the algorithm
 *
 * @author Martin Guetlein, Andreas Maunz
 * @version 1.0 20/9/2012
 */
public class GetCDKDescriptors {
	private static DescriptorEngine ENGINE = new DescriptorEngine(DescriptorEngine.ATOMIC);


  /**
  * Example main
  *
  */
  public static void main(String args[]) throws java.io.IOException 
  {
    String inpath = "HMDB00001.sdf";
    String outpath = "HMDB00001.csv";
    getDescriptorCSV(inpath,outpath,"");
  }
 /**
  * Constructor, executing the algorithm
  *
  * @param: string The path to the input SDF file
  * @param: string The path to the output CSV file
  */
  public GetCDKDescriptors(String inpath, String outpath, String descNamesStr) throws java.io.IOException  { 
    getDescriptorCSV(inpath,outpath,descNamesStr);
  }



 /**
  * Calculate descriptors. Omits IPMolecularLearningDescriptor
  *
  * @param string path to SDF input file
  * @param string path to CSV output file
  * @param string comma-seperated list of descriptor names (if empty, all descriptors will be calculated)
  */
 public static void getDescriptorCSV(String sdfInputPath, String csvOutputPath, String descNamesStr) throws java.io.IOException  
 {
    List<IMolecule> mols = readMolecules(sdfInputPath);
		System.err.println("read " + mols.size() + " compounds");
		List<IDescriptor> descriptors = ENGINE.getDescriptorInstances();
		System.err.println("found " + descriptors.size() + " descriptors");

    List<String> descNames = Arrays.asList(descNamesStr.split(","));
    ArrayList<String> colNames = new ArrayList<String>();
    ArrayList<Double[]> values = new ArrayList<Double[]>();
    for (IDescriptor desc : descriptors) {
      if (desc instanceof IPAtomicLearningDescriptor)
        continue;
      String tname = desc.getClass().getName();
      String[] tnamebits = tname.split("\\.");
      tname = tnamebits[tnamebits.length-1];
      if ((descNamesStr.length()>0) && (!descNames.contains(tname)))
        continue;
      String[] colNamesArr = desc.getDescriptorNames();
      for (int idx=0; idx<colNamesArr.length; idx++) {
        colNamesArr[idx] = tname + "-" + colNamesArr[idx];
      }
      colNames.addAll(Arrays.asList(colNamesArr));
      for (IAtomContainer mol : mols) {
        int atomCount = mol.getAtomCount();
        List<IAtom> atoms = new ArrayList<IAtom>();
        for (int i = 0; i < atomCount; i++) {
          atoms.add(mol.getAtom(i));
        }
        values.addAll(computeListsAtomic(mol, atoms, (IAtomicDescriptor) desc));
      }
    }

    int ncol = values.size();
    int nrow = mols.get(0).getAtomCount();
    FileWriter fstream = new FileWriter(csvOutputPath);
    BufferedWriter out = new BufferedWriter(fstream);
    out.write("SMILES,");
    for (int c=0; c<ncol; c++) {
      if (c!=0) out.write(",");
      out.write(colNames.get(c));
    }
    out.write("\n");
    for (int r=0; r<nrow; r++) {
      IAtom atom = mols.get(0).getAtom(r);
      String smi = atom.getSymbol();
      out.write(smi + ",");
      for (int c=0; c<ncol; c++) {
        if (c!=0) out.write(",");
        out.write(""+values.get(c)[r]);
      }
      out.write("\n");
    }
    out.flush();
 }

 public static ArrayList<Double[]> getAtomicDescriptor(String sdf, String descNamesStr) throws java.io.IOException  
 {
    List<IMolecule> mols = readMoleculesString(sdf);
    System.err.println("read " + mols.size() + " compounds");
    List<IDescriptor> descriptors = ENGINE.getDescriptorInstances();
    System.err.println("found " + descriptors.size() + " descriptors");

    List<String> descNames = Arrays.asList(descNamesStr.split(","));
    ArrayList<String> colNames = new ArrayList<String>();
    ArrayList<Double[]> values = new ArrayList<Double[]>();
    for (IDescriptor desc : descriptors) {
      if (desc instanceof IPAtomicLearningDescriptor)
        continue;
      String tname = desc.getClass().getName();
      String[] tnamebits = tname.split("\\.");
      tname = tnamebits[tnamebits.length-1];
      if ((descNamesStr.length()>0) && (!descNames.contains(tname)))
        continue;
      String[] colNamesArr = desc.getDescriptorNames();
      for (int idx=0; idx<colNamesArr.length; idx++) {
        colNamesArr[idx] = tname + "-" + colNamesArr[idx];
      }
      colNames.addAll(Arrays.asList(colNamesArr));
      for (IAtomContainer mol : mols) {
        int atomCount = mol.getAtomCount();
        List<IAtom> atoms = new ArrayList<IAtom>();
        for (int i = 0; i < atomCount; i++) {
          atoms.add(mol.getAtom(i));
        }
        values.addAll(computeListsAtomic(mol, atoms, (IAtomicDescriptor) desc));

        // for (int i = 0; i < atomCount; i++) {
        //   System.out.println(mol.getAtom(i).getSymbol());
        // }
      }
    }

    return values;
 }

 // Find all relevant hydrogen atoms in molecule
public static ArrayList<String> getHydrogenAtoms(String sdf) {
  List<IMolecule> mols = readMoleculesString(sdf);
  ArrayList<String> hydrogen_positions = new ArrayList<String>();

  for (IAtomContainer mol : mols) {
    int atomCount = mol.getAtomCount();
    for (int i = 0; i < atomCount; i++) {
      IAtom atom = mol.getAtom(i);

      if (atom.getSymbol().equals("H")) {
        List<IAtom> connected_atoms = mol.getConnectedAtomsList(atom);
        if (!connected_atoms.get(0).getSymbol().equals("O") && !connected_atoms.get(0).getSymbol().equals("N")) {
          hydrogen_positions.add(String.valueOf(i));
        }
      }
    }
  }
  return hydrogen_positions;
}

// Find nearest atom to all atoms in a molecule
public static ArrayList<ArrayList<String>> getNearestAtoms(String sdf) {
  List<IMolecule> mols = readMoleculesString(sdf);
  ArrayList<ArrayList<String>> atom_distances = new ArrayList<ArrayList<String>>();

  for (IAtomContainer mol : mols) {
    int atomCount = mol.getAtomCount();
    List<IAtom> atoms = new ArrayList<IAtom>();
    for (int i = 0; i < atomCount; i++) {
      atoms.add(mol.getAtom(i));
    }
    for (int i = 0; i < atoms.size(); i++)
    {
      Double[] distances = new Double[atoms.size()];
      // Double min_distance = 999999.9999;
      // int min_d_index = -1;
      for (int j = 0; j < atoms.size(); j++) {
        if (j == i) {
          // Large number so that sorting puts it last
          distances[j] = 99999.12;
          continue;
        }

        Point3d firstPoint = atoms.get(i).getPoint3d();
        Point3d secondPoint = atoms.get(j).getPoint3d();
        Double distance = firstPoint.distance(secondPoint);
        distances[j] = distance;

        // if (distance < min_distance) {
        //   min_distance = distance;
        //   min_d_index = j;
        // }
      }

      //atom_distances.add(String.valueOf(min_d_index));

      ArrayList<String> indices = new ArrayList<String>();
      Double[] d = distances.clone();
      Arrays.sort(d);
      List<Double> d_list = Arrays.asList(distances);

      for (int j = 0; j < distances.length; j++)
      {
        String index = String.valueOf(d_list.indexOf(d[j]));
        indices.add(index);
      }
      atom_distances.add(indices);
    }

  }
  return atom_distances;
}

// Calculate RDF Proton descriptors by calculating aromaticity

// public static ArrayList<String> calculateRDFProtonDescriptor(IAtomContainer molecule, Descriptor desc) {
//   ElectronDonation model       = ElectronDonation.cdk();
//   CycleFinder      cycles      = Cycles.cdkAromaticSet();
//   Aromaticity      aromaticity = new Aromaticity(model, cycles);

//   // apply our configured model to each molecule, the CDK model
//   // requires that atom types are perceived
//   AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
//   boolean aromaticity = aromaticity.apply(molecule);

//   int atomCount = mol.getAtomCount();
//   List<IAtom> atoms = new ArrayList<IAtom>();

//   for (int i = 0; i < atomCount; i++) {
//     atoms.add(mol.getAtom(i));
//   }
//   // check for RDF descriptors
//   Array<String> value_desc = new Array<String>();
//   for (int i = 0; i < atoms.size(); i++) {
//     value_desc.add(Double.valueOf(descriptor.calculate(molecule, atoms.get(i), aromaticity));
//   }
// }


 /**
  * Get SMILES code for a molecule
  *
  * @param IMolecule The molecule
  * @return string The SMILES code
  */
 public static String getSmiles(IMolecule m)
 {
   Map<Object, Object> props = m.getProperties();
   for (Object key : props.keySet()) {
     if (key.toString().equals("STRUCTURE_SMILES") || key.toString().equals("SMILES"))
       return props.get(key).toString();
   }
   SmilesGenerator g = new SmilesGenerator();
   return g.createSMILES(m);
 }


	public static List<Double[]> computeLists(List<IMolecule> mols, IMolecularDescriptor desc )
	{
    System.out.println("computing descriptor " + getName(desc));
    List<Double[]> values = computeDescriptors(mols, (IMolecularDescriptor) desc);
    return values;
	}

  public static List<Double[]> computeListsAtomic(IAtomContainer mol, List<IAtom> atoms, IAtomicDescriptor desc )
  {
    System.out.println("computing descriptor " + getName(desc));
    List<Double[]> values = computeDescriptorsAtomic(mol, atoms, desc);
    return values;
  }


 /**
  * Read in molecules, using any supported format
  *
  * @param string The input file
  * @return Vector<IMolecule> The molecules
  */
	public static List<IMolecule> readMolecules(String filepath)
	{
		Vector<IMolecule> mols = new Vector<IMolecule>();
		File file = new File(filepath);
		if (!file.exists())
			throw new IllegalArgumentException("file not found: " + filepath);
		List<IAtomContainer> list;
		try
		{
			ISimpleChemObjectReader reader = new ReaderFactory().createReader(new InputStreamReader(
					new FileInputStream(file)));
			if (reader == null)
				throw new IllegalArgumentException("Could not determine input file type");
			IChemFile content = (IChemFile) reader.read((IChemObject) new ChemFile());
			list = ChemFileManipulator.getAllAtomContainers(content);
			reader.close();
		}
		catch (Exception e)
		{
			e.printStackTrace();
			return null;
		}

		for (IAtomContainer iAtomContainer : list)
		{
			IMolecule mol = (IMolecule) iAtomContainer;
			mol = (IMolecule) AtomContainerManipulator.removeHydrogens(mol);
			try
			{
				AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
			}
			catch (Exception e)
			{
				e.printStackTrace();
			}
			try
			{
				CDKHueckelAromaticityDetector.detectAromaticity(mol);
			}
			catch (Exception e)
			{
				e.printStackTrace();
			}
			if (mol.getAtomCount() == 0)
				System.err.println("molecule has no atoms");
			else
				mols.add(mol);
		}
		return mols;
	}

  public static List<IMolecule> readMoleculesString(String sdf)
  {
    Vector<IMolecule> mols = new Vector<IMolecule>();
    if (sdf.equals(""))
      throw new IllegalArgumentException("No sdf found" + sdf);
    List<IAtomContainer> list;
    try
    {
      InputStream is = new ByteArrayInputStream(sdf.getBytes("UTF-8"));
      ISimpleChemObjectReader reader = new ReaderFactory().createReader(new InputStreamReader(is));
      if (reader == null)
        throw new IllegalArgumentException("Could not determine input file type");
      IChemFile content = (IChemFile) reader.read((IChemObject) new ChemFile());
      list = ChemFileManipulator.getAllAtomContainers(content);
      reader.close();
    }
    catch (Exception e)
    {
      e.printStackTrace();
      return null;
    }

    for (IAtomContainer iAtomContainer : list)
    {
      IMolecule mol = (IMolecule) iAtomContainer;
      //mol = (IMolecule) AtomContainerManipulator.removeHydrogens(mol);
      try
      {
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
      }
      catch (Exception e)
      {
        e.printStackTrace();
      }
      try
      {
        CDKHueckelAromaticityDetector.detectAromaticity(mol);
      }
      catch (Exception e)
      {
        e.printStackTrace();
      }
      if (mol.getAtomCount() == 0)
        System.err.println("molecule has no atoms");
      else
        mols.add(mol);
    }
    return mols;
  }

	public static List<Double[]> computeDescriptors(List<IMolecule> mols, IMolecularDescriptor descriptor)
	{
		List<Double[]> vv = new ArrayList<Double[]>();

		for (int j = 0; j < getSize(descriptor); j++)
			vv.add(new Double[mols.size()]);

		for (int i = 0; i < mols.size(); i++)
		{
			if (mols.get(i).getAtomCount() == 0)
			{
				for (int j = 0; j < getSize(descriptor); j++)
					vv.get(j)[i] = null;
			}
			else
			{
				try
				{
					IDescriptorResult res = descriptor.calculate(mols.get(i)).getValue();
					if (res instanceof IntegerResult)
						vv.get(0)[i] = (double) ((IntegerResult) res).intValue();
					else if (res instanceof DoubleResult)
						vv.get(0)[i] = ((DoubleResult) res).doubleValue();
					else if (res instanceof DoubleArrayResult)
						for (int j = 0; j < getSize(descriptor); j++)
							vv.get(j)[i] = ((DoubleArrayResult) res).get(j);
					else if (res instanceof IntegerArrayResult)
						for (int j = 0; j < getSize(descriptor); j++)
							vv.get(j)[i] = (double) ((IntegerArrayResult) res).get(j);
					else
						throw new IllegalStateException("Unknown idescriptor result value for '" + descriptor + "' : "
								+ res.getClass());
				}
				catch (Throwable e)
				{
					System.err.println("Could not compute cdk feature " + descriptor);
					e.printStackTrace();
					for (int j = 0; j < getSize(descriptor); j++)
						vv.get(j)[i] = null;
				}
			}
			for (int j = 0; j < getSize(descriptor); j++)
				if (vv.get(j)[i] != null && (vv.get(j)[i].isNaN() || vv.get(j)[i].isInfinite()))
					vv.get(j)[i] = null;
		}

		return vv;
	}

  public static List<Double[]> computeDescriptorsAtomic(IAtomContainer mol, List<IAtom> atoms, IAtomicDescriptor descriptor)
  {
    List<Double[]> vv = new ArrayList<Double[]>();
    vv.add(new Double[atoms.size()]);
    for (int i = 0; i < atoms.size(); i++)
    {
      if (atoms.get(i) == null)
      {
          vv.get(0)[i] = null;
      }
      else
      {
        try
        {
          IDescriptorResult res = descriptor.calculate(atoms.get(i), mol).getValue();
          if (res instanceof IntegerResult) {
            vv.get(0)[i] = (double) ((IntegerResult) res).intValue();
            //System.out.println(vv.get(0)[i]); 
          }
          else if (res instanceof DoubleResult) {
            vv.get(0)[i] = ((DoubleResult) res).doubleValue();
            //System.out.println(vv.get(0)[i]); 
          }
          else if (res instanceof DoubleArrayResult) {
            vv.get(0)[i] = ((DoubleArrayResult) res).get(0); 
            //System.out.println(vv.get(0)[i]);  
          }
          else if (res instanceof IntegerArrayResult) {
            vv.get(0)[i] = (double) ((IntegerArrayResult) res).get(0);
            //System.out.println(vv.get(0)[i]); 
          }
          else
            throw new IllegalStateException("Unknown idescriptor result value for '" + descriptor + "' : "
                + res.getClass());
        }
        catch (Throwable e)
        {
          System.err.println("Could not compute cdk feature " + descriptor);
          e.printStackTrace();
            vv.get(0)[i] = 0.0;
        }
      }
      if (vv.get(0)[i] != null && (vv.get(0)[i].isNaN() || vv.get(0)[i].isInfinite()))
        vv.get(0)[i] = 0.0;
    }

    return vv;
  }


 /**
  * Get length of result for a given descriptor
  *
  * @param IMolecularDescriptor The descriptor
  * @return int The length
  */
	private static int getSize(IMolecularDescriptor descriptor) 
  {
		IDescriptorResult r = descriptor.getDescriptorResultType();
		if (r instanceof DoubleArrayResultType)
			return ((DoubleArrayResultType) r).length();
		else if (r instanceof IntegerArrayResultType)
			return ((IntegerArrayResultType) r).length();
		else
			return 1;
	}


 /**
  * Get name for a given descriptor
  *
  * @param IMolecularDescriptor The descriptor
  */
	private static String getName(IDescriptor descriptor) 
  {
    try
    {
      String name = ENGINE.getDictionaryTitle(descriptor.getSpecification()).trim();
  		if (name != null) {
        return name;
      }
      else {
        return "";
      }
    }
    catch (Throwable e)
    {
      return "";
    }
	}


}
