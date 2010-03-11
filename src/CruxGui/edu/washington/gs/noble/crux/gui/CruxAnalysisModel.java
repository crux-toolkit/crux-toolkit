package edu.washington.gs.noble.crux.gui;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.io.PrintWriter;
import java.io.Serializable;
import java.io.Writer;
import java.util.logging.Logger;

import javax.swing.JOptionPane;

public class CruxAnalysisModel extends Object implements Serializable{
	
	private static final long serialVersionUID = 3352659336939985777L;
	private static Logger logger = 
		Logger.getLogger("edu.washington.gs.noble.crux.gui");


	public enum CruxComponents {
		CREATE_INDEX("create-index"), 
		SEARCH_FOR_MATCHES("search-for-matches");
		String name;
		CruxComponents(String name) {this.name = name;}
		public String toString() { return this.name; }
		public static int getSize() { return 2; };
	};
	
	public enum IsotopicMassType { 
		AVERAGE("average"), 
		MONOISOTOPIC("mono");
		String name;
		IsotopicMassType(String name) {this.name = name;}
		public String toString() { return this.name; }
		public static int getSize() { return 2; };
	};
	
	public enum DigestType {
		FULL("full-digest"), 
		PARTIAL("partial-digest"); 
		String name;
		DigestType(String name) {this.name = name;}
		public String toString() { return this.name; }
		public static int getSize() { return 2; };
	};
	
	public enum Enzyme {
		TRYPSIN("trypsin"),
		CHYMOTRYPSIN("chymotrypsin"), 
		ELASTASE("elastase"), 
		CLOSTRIPAIN("clostripain"),
		CYANOGEN("cyanogen-bromide"),
		IDOSOBENZOATE("idosobenzoate"),
		PROLINE_ENDO("proline-endopeptidase"),
		STAPH_PROT("staph-protease"),
		MOD_CHYMOTRYPSIN("modified-chymotrypsin"),
		ETC("elastase-trypsin-chymotrypsin"),
		NONE("no-enzyme");
		String name;
		Enzyme(String name) {this.name = name;}
		public static int getSize() { return 11; };
		public String toString() { return this.name; };
	}
	
	private String name;
	private boolean componentsToRun[]= new boolean[CruxComponents.getSize()];
	private String proteinSource;
	private String spectraSource;
	private IsotopicMassType massType;
	private DigestType digestType;
	private Enzyme enzyme;
	private boolean allowMissedCleavages;
	
	public CruxAnalysisModel() {
		restoreDefaults();
	}
	
	String getName() {
		return name;
	}
	
	void setName(String name) {
		this.name = name;
		logger.info("Model parameter 'name' set to " + name);
	}
	
	void restoreDefaults() {
		name = null;
		componentsToRun[CruxComponents.CREATE_INDEX.ordinal()] = false;
		componentsToRun[CruxComponents.SEARCH_FOR_MATCHES.ordinal()] = false;
		proteinSource = null;
		spectraSource = null;
		massType = IsotopicMassType.AVERAGE;
		digestType = DigestType.FULL;
		enzyme = Enzyme.TRYPSIN;
		allowMissedCleavages = false;
		logger.info("Model restored to default values");
	}
	
	boolean getRunCreateIndex() {
		return componentsToRun[CruxComponents.CREATE_INDEX.ordinal()];
	}
	
	void setRunCreateIndex(boolean value) {
		componentsToRun[CruxComponents.CREATE_INDEX.ordinal()] = value;
		logger.info("Run 'create-index' set to " + value);
	}
	
	boolean getRunSearchForMatches() {
		return componentsToRun[CruxComponents.SEARCH_FOR_MATCHES.ordinal()];
	}

	void setRunSearchForMatches(boolean value) {
		componentsToRun[CruxComponents.SEARCH_FOR_MATCHES.ordinal()] = value;
		logger.info("Run 'search-for-matches' set to " + value);
	}
	
	IsotopicMassType getMassType() {
		return this.massType;
	}

	void setMassType(IsotopicMassType type) {
		this.massType = type;
		logger.info("Model parameter 'massType' set to " + type.toString());
	}
	
	DigestType getDigestType() {
		return this.digestType;
	}
	
	void setDigestType(DigestType type) {
		this.digestType =type;
		logger.info("Model parameter 'digestType' set to " + type.toString());
	}
	
	Enzyme getEnzyme() {
		return this.enzyme;
	}

	void setEnzyme(Enzyme e) {
		this.enzyme = e;
		logger.info("Model parameter 'enzyme' set to " + e.toString());
	}
	
	boolean getAllowMissedCleavages() {
		return this.allowMissedCleavages;
	}
	
	void setAllowMissedCleavages(boolean value) {
		this.allowMissedCleavages = value;
		logger.info("Model parameter 'allowMissedCleavages' set to " + value);
	}

	void setProteinDatabase(String database) {
		this.proteinSource = database;
		logger.info("Model parameter 'proteinDatabase' set to " + database);
	}
	
	public boolean toBinaryFile() {
		// Serialze the model to the analysis directory
		String model_path = name + "/" + name + ".model";
		boolean result = false;
		try {
			FileOutputStream fileOut = new FileOutputStream(model_path);
			ObjectOutputStream out = new ObjectOutputStream(fileOut);
			out.writeObject(this);
			out.close();
			logger.info("Saved model to " + model_path);	
			result = true;
		} catch (IOException e) {
			logger.info("Unable to saved model: " + e.toString());
			JOptionPane.showMessageDialog(null, "Unable to save model: " + e.toString());
		}
		return result;
	}

	public boolean toParameterFile() {
		boolean result = false;
		String fileName = name + "/" + "analysis.params";
		try {
			PrintWriter out = new PrintWriter(new File(fileName));
			out.println("#Parameter file for analysis " + name);
			out.println("isotopic-mass=" + massType);
			out.println("digestion=" + digestType);
			out.println("enzyme=" + enzyme);
			if (allowMissedCleavages) {
				out.println("missed-cleavages=T");
			}
			else {
				out.println("missed-cleavages=F");
			}
			out.close();
			logger.info("Saved analysis parameter file to " + fileName);	
			result = true;
		}
		catch (IOException e) {
			logger.info("Unable to write parameter file for analsyis: " + e.toString());
			JOptionPane.showMessageDialog(null, "Unable to write parameter file for analsyis: " + e.toString());
		}
		return result;
	}

	boolean isValidAnalysis() {
		
		boolean result = true;
		
		if (name == null) {
			JOptionPane.showMessageDialog(null, "Please choose a name for this analysis.");
			result = false;
		}
		if (componentsToRun[CruxComponents.CREATE_INDEX.ordinal()]) {
			// Must have specified protein database
			if (proteinSource == null) {
				JOptionPane.showMessageDialog(null, "A protein database must be specified to create an index.");
				result = false;
			}
		}
		if (componentsToRun[CruxComponents.SEARCH_FOR_MATCHES.ordinal()]) {
			// Must have specified protein database and spectra source
			if (proteinSource == null || spectraSource == null) {
				JOptionPane.showMessageDialog(null, "A protein database and a spectra source must be specified to search for maatches.");
				result = false;
			}
		}
		logger.info("Validation of analysis: " + result + ".");

		return result;
	}
	
	boolean run() {
		
		if (!isValidAnalysis()) {
			return false;
		}
		
		// Write out binary file of the analysis model to the analysis directory
		if (!toBinaryFile()) {
			return false;
		}
		
		// Write out a parameter file to the analysis directory
		if (!toParameterFile()) {
			return false;
		}
	
		if (componentsToRun[CruxComponents.CREATE_INDEX.ordinal()]) {
			//  Run the component using the component directory and parameter file
			String command = "crux create-index --parameter-file " + name + "/analysis.params " + proteinSource + " " + name + "/" + CruxComponents.CREATE_INDEX.toString();
			try {
				logger.info("Attempte to execute command: " + command);
				Process process = Runtime.getRuntime().exec(command);
			}
			catch (IOException e) {
				logger.info("Unable to execute command: " + e.toString());
				JOptionPane.showMessageDialog(null, "Unable to write parameter file for analsyis: " + e.toString());
				return false;
			}
		}
		if (componentsToRun[CruxComponents.SEARCH_FOR_MATCHES.ordinal()]) {
			//	create a subdirectory in the analysis directory
			//  write out a parameter file in the component directory
			//  run the component using the component directory and parameter file
		}

		return true;
	}
}