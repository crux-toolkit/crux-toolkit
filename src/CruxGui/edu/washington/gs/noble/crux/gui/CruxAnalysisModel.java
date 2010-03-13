package edu.washington.gs.noble.crux.gui;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.io.PrintWriter;
import java.io.Serializable;
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
	private int numComponentsToRun;
	private String proteinSource;
	private String spectraSource;
	private IsotopicMassType massType;
	private DigestType digestType;
	private Enzyme enzyme;
	private boolean allowMissedCleavages;
	
	public CruxAnalysisModel() {
		restoreDefaults();
	}
	
	public String getName() {
		return name;
	}
	
	public void setName(String name) {
		this.name = name;
		logger.info("Model parameter 'name' set to " + name);
	}
	
	public void restoreDefaults() {
		name = null;
		numComponentsToRun = 0;
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
	
	public boolean getRunCreateIndex() {
		return componentsToRun[CruxComponents.CREATE_INDEX.ordinal()];
	}
	
	public void setRunCreateIndex(boolean value) {
		componentsToRun[CruxComponents.CREATE_INDEX.ordinal()] = value;
		if (value) {
		    ++numComponentsToRun;
		}
		else {
		    --numComponentsToRun;
		}
		logger.info("Run 'create-index' set to " + value);
	}
	
	public boolean getRunSearchForMatches() {
		return componentsToRun[CruxComponents.SEARCH_FOR_MATCHES.ordinal()];
	}

	public void setRunSearchForMatches(boolean value) {
		componentsToRun[CruxComponents.SEARCH_FOR_MATCHES.ordinal()] = value;
		if (value) {
		    ++numComponentsToRun;
		}
		else {
		    --numComponentsToRun;
		}
		logger.info("Run 'search-for-matches' set to " + value);
	}
	
	public IsotopicMassType getMassType() {
		return this.massType;
	}

	public void setMassType(IsotopicMassType type) {
		this.massType = type;
		logger.info("Model parameter 'massType' set to " + type.toString());
	}
	
	public DigestType getDigestType() {
		return this.digestType;
	}
	
	public void setDigestType(DigestType type) {
		this.digestType =type;
		logger.info("Model parameter 'digestType' set to " + type.toString());
	}
	
	public Enzyme getEnzyme() {
		return this.enzyme;
	}

	public void setEnzyme(Enzyme e) {
		this.enzyme = e;
		logger.info("Model parameter 'enzyme' set to " + e.toString());
	}
	
	public boolean getAllowMissedCleavages() {
		return this.allowMissedCleavages;
	}
	
	public void setAllowMissedCleavages(boolean value) {
		this.allowMissedCleavages = value;
		logger.info("Model parameter 'allowMissedCleavages' set to " + value);
	}

	public String getProteinDatabase() {
		return proteinSource;
	}
	
	public void setProteinDatabase(String database) {
		this.proteinSource = database;
		logger.info("Model parameter 'proteinDatabase' set to " + database);
	}
	
	public String getSpectraFilemame() {
		return spectraSource;
	}
	
	public void setSpectraFilename(String fileName) {
		this.spectraSource = fileName;
		logger.info("Model parameter 'spectraSource' set to " + fileName);
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
		
		if (numComponentsToRun <= 0) {
			JOptionPane.showMessageDialog(null, "Please select some components to run.");
			result = false;
		}
		if (name == null) {
			JOptionPane.showMessageDialog(null, "Please choose a name for this analysis.");
			result = false;
		}
		if (componentsToRun[CruxComponents.CREATE_INDEX.ordinal()]) {
			// Must have specified protein database
			if (proteinSource == null || proteinSource.length() == 0) {
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
	
	Process run() {
		
		Process process = null;
		
		if (!isValidAnalysis()) {
			return process;
		}
		
		// Write out binary file of the analysis model to the analysis directory
		if (!toBinaryFile()) {
			return process;
		}
		
		// Write out a parameter file to the analysis directory
		if (!toParameterFile()) {
			return process;
		}
	
		if (componentsToRun[CruxComponents.CREATE_INDEX.ordinal()]) {
			//  Run the component using the component directory and parameter file
			String subCommand = CruxComponents.CREATE_INDEX.toString();
			String command = "crux " + subCommand + " --parameter-file " + name + "/analysis.params " + proteinSource + " " + name + "/" + subCommand;
			try {
				logger.info("Attempte to execute command: " + command);
				process = Runtime.getRuntime().exec(command);
			}
			catch (IOException e) {
				logger.info("Unable to execute command: " + e.toString());
				JOptionPane.showMessageDialog(null, "Unable to execute command: " + e.toString());
				return process;
			}
		}
		if (componentsToRun[CruxComponents.SEARCH_FOR_MATCHES.ordinal()]) {
			//  Run the component using the component directory and parameter file
			String proteinSource;
			if (componentsToRun[CruxComponents.CREATE_INDEX.ordinal()]) {
				proteinSource = name + "/" + CruxComponents.CREATE_INDEX.toString();
			}
			else {
				proteinSource = this.proteinSource;
			}
			String subCommand = CruxComponents.SEARCH_FOR_MATCHES.toString();
			String command = "crux " + subCommand + " --overwrite T --output-dir " + name + "/" + subCommand + " --parameter-file " + name + "/analysis.params " + spectraSource + " " + proteinSource;
			try {
				logger.info("Attempte to execute command: " + command);
				process = Runtime.getRuntime().exec(command);
			}
			catch (IOException e) {
				logger.info("Unable to execute command: " + e.toString());
				JOptionPane.showMessageDialog(null, "Unable to execute command: " + e.toString());
				return process;
			}
		}

		return process;
	}
}