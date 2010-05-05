package edu.washington.gs.noble.crux.gui;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.logging.Logger;

import javax.swing.JOptionPane;

/**
 * The CruxAnalysisClass is used to hold all the information required to run a Crux analysis.
 * This includes which crux components to run, the status of the run, the parameters for those 
 * tools, the location of the crux binary, the name of the analysis, and names of the directories 
 * containing the output of the analysis. Instances of the class can be serialized to the file system
 * to provide permanent storage for the details of the analysis.
 * 
 * @author Charles E. Grant
 *
 */

public class CruxAnalysisModel extends Object implements Serializable{
	
	private static final long serialVersionUID = 3352659336939985777L;
	private static Logger logger = 
		Logger.getLogger("edu.washington.gs.noble.crux.gui");


	/** The list of the component tools of crux */
	public enum CruxComponents {
		CREATE_INDEX("create-index"), 
		SEARCH_FOR_MATCHES("search-for-matches"),
		COMPUTE_Q_VALUES("compute-q-values"),
		PERCOLATOR("percolator"),
		QRANKER("qranker");
		String name;
		CruxComponents(String name) {this.name = name;}
		public String toString() { return this.name; }
		public static int getSize() { return 5; };
	};
	
	/** The list of run states */
	public enum RunStatus {
		NOT_RUN("not-run"), 
		COMPLETED("completed"), 
		FAILED("failed");
		String name;
		RunStatus(String name) {this.name = name;}
		public String toString() { return this.name; }
		public static int getSize() { return 3; };
	};

	/** The list of allowed isotopic mass types */
	public enum IsotopicMassType { 
		AVERAGE("average"), 
		MONOISOTOPIC("mono");
		String name;
		IsotopicMassType(String name) {this.name = name;}
		public String toString() { return this.name; }
		public static int getSize() { return 2; };
	};
	
	/** The list of allowed peptide digestion types */
	public enum DigestType {
		FULL("full-digest"), 
		PARTIAL("partial-digest"); 
		String name;
		DigestType(String name) {this.name = name;}
		public String toString() { return this.name; }
		public static int getSize() { return 2; };
	};
	
	/** The list of allowed enzymes */
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
	
	/** The list of allowed spectrum charges */
	public enum AllowedSpectrumCharge {
		ONE("1"), 
		TWO("2"),
		THREE("3"),
		ALL("all") ; 
		String name;
		AllowedSpectrumCharge(String name) {this.name = name;}
		public String toString() { return this.name; }
		public static int getSize() { return 4; };
	};
	
	private boolean needsSaving = false;
	private String name;
	private String pathToCrux;
	private final boolean componentsToRun[]= new boolean[CruxComponents.getSize()];
	private final RunStatus componentsRunStatus[]= new RunStatus[CruxComponents.getSize()];
	private final boolean showAdvancedParameters[]= new boolean[CruxComponents.getSize()];
	private int numComponentsToRun;
	
	// Crux parameters
	private boolean allowMissedCleavages;
	private String customEnzymeBeforeCleavage;
	private String customEnzymeAfterCleavage;
	private boolean decoyPValues;
	private DigestType digestType;
	private Enzyme enzyme;
	private boolean featureFile;
	private IsotopicMassType massType;
	private int maxLength;
	private double maxMass;
	private int maxMods;
	private int minLength;
	private double minMass;
	private int numDecoysPerTarget;
	private int numDecoyFiles;
	private double pi0;
	private int printSearchProgress;
	private String proteinSource;
	private double seed;
	private String spectraSource;
	private AllowedSpectrumCharge spectrumCharge;
	private double spectrumMaxMass;
	private double spectrumMinMass;
	private int topMatch;
	
	// Crux parameters defaults
	final private boolean allowMissedCleavagesDefault = false;
	final private String customEnzymeAfterCleavageDefault = null;
	final private String customEnzymeBeforeCleavageDefault = null;
	final private boolean decoyPValuesDefault = false;
	final private DigestType digestTypeDefault = DigestType.FULL;
	final private Enzyme enzymeDefault = Enzyme.TRYPSIN;
	final private boolean featureFileDefault = false;
	final private IsotopicMassType massTypeDefault = IsotopicMassType.AVERAGE;
	final private int maxLengthDefault = 50;
	final private double maxMassDefault = 7200;
	final private int maxModsDefaults = -1;
	final private int minLengthDefault = 6;
	final private double minMassDefault = 200;
	final private int numDecoyFilesDefault = 2;
	final private int numDecoysPerTargetDefault = 2;
	final private double pi0Default = 0.9;
	final private int printSearchProgressDefault = 10;
	final private String proteinSourceDefault = null;
	final private String spectraSourceDefault = null;
	final private AllowedSpectrumCharge spectrumChargeDefault = AllowedSpectrumCharge.ALL;
	final private double seedDefault = -1.0;
	final private double spectrumMaxMassDefault = -1.0;
	final private double spectrumMinMassDefault = 0.0;
	final private int topMatchDefault = 5;
	
	public CruxAnalysisModel() {
		setDefaults();
	}
	
	/** Serialize the analysis to the analysis directory */
	public boolean saveModelToBinaryFile() {
		String model_path = name + "/" + name + ".model";
		boolean result = false;
		boolean oldNeedsSaving = needsSaving();
		try {
			setNeedsSaving(false);
			FileOutputStream fileOut = new FileOutputStream(model_path);
			ObjectOutputStream out = new ObjectOutputStream(fileOut);
			out.writeObject(this);
			out.close();
			logger.info("Saved model to " + model_path);	
			result = true;
		} catch (IOException e) {
			logger.info("Unable to saved mode: " + e.toString());
			JOptionPane.showMessageDialog(null, "Unable to save model: " + e.toString());
		}
		return result;
	}
	
	public CruxAnalysisModel readModelFromBinaryFile(String name) {
		CruxAnalysisModel model = null;
		String model_path = name + "/" + name + ".model";
		try {
			FileInputStream fileIn = new FileInputStream(model_path);
			ObjectInputStream in = new ObjectInputStream(fileIn);
			model = (CruxAnalysisModel) in.readObject();
		} catch (IOException e) {
			logger.info("Unable to load model " + name + ": " + e.toString());
			JOptionPane.showMessageDialog(null, "Unable to load model " + name + ": "  + e.toString());
		} catch (ClassNotFoundException e) {
			logger.info("Unable to load model " + name + ": " + e.toString());
			JOptionPane.showMessageDialog(null, "Unable to load model " + name + ": "  + e.toString());
		}
		return model;
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

	boolean needsSaving() {
		return needsSaving;
	}
	
	void setNeedsSaving(boolean needsSaving) {
	    this.needsSaving = needsSaving;
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
		if (needsSaving) {
    		if (!saveModelToBinaryFile()) {
    			return process;
    		}
		}
		
		// Write out a parameter file to the analysis directory
		if (!toParameterFile()) {
			return process;
		}
	
		if (componentsToRun[CruxComponents.CREATE_INDEX.ordinal()] 
		    && componentsRunStatus[CruxComponents.CREATE_INDEX.ordinal()] == RunStatus.NOT_RUN) {
			//  Run the component using the component directory and parameter file
			String subCommand = CruxComponents.CREATE_INDEX.toString();
			String command = pathToCrux + " " + subCommand + " --parameter-file " + name + "/analysis.params " + proteinSource + " " + name + "/" + subCommand;
			try {
				logger.info("Attempte to execute command: " + command);
				process = Runtime.getRuntime().exec(command);
				process.waitFor();
				if (process.exitValue() == 0) {
				    componentsRunStatus[CruxComponents.CREATE_INDEX.ordinal()] = RunStatus.COMPLETED;
				}
				else {
				    componentsRunStatus[CruxComponents.CREATE_INDEX.ordinal()] = RunStatus.FAILED;
				}
			} catch (Exception e) {
				logger.info("Unable to execute command: " + e.toString());
				JOptionPane.showMessageDialog(null, "Unable to execute command: " + e.toString());
				return process;
			}
		}
		if (componentsToRun[CruxComponents.SEARCH_FOR_MATCHES.ordinal()]
		    && componentsRunStatus[CruxComponents.SEARCH_FOR_MATCHES.ordinal()] == RunStatus.NOT_RUN) {
			//  Run the component using the component directory and parameter file
			String proteinSource;
			if (componentsToRun[CruxComponents.CREATE_INDEX.ordinal()]) {
				proteinSource = name + "/" + CruxComponents.CREATE_INDEX.toString();
			}
			else {
				proteinSource = this.proteinSource;
			}
			String subCommand = CruxComponents.SEARCH_FOR_MATCHES.toString();
			String command = pathToCrux + " " + subCommand + " --overwrite T --output-dir " + name + "/" + subCommand + " --parameter-file " + name + "/analysis.params " + spectraSource + " " + proteinSource;
			try {
				logger.info("Attempte to execute command: " + command);
				process = Runtime.getRuntime().exec(command);
				process.waitFor();
				if (process.exitValue() == 0) {
				    componentsRunStatus[CruxComponents.SEARCH_FOR_MATCHES.ordinal()] = RunStatus.COMPLETED;
				}
				else {
				    componentsRunStatus[CruxComponents.SEARCH_FOR_MATCHES.ordinal()] = RunStatus.FAILED;
				}
			}
			catch (Exception e) {
				logger.info("Unable to execute command: " + e.toString());
				JOptionPane.showMessageDialog(null, "Unable to execute command: " + e.toString());
				return process;
			}
		}

		return process;
	}

	public String getName() {
		return name;
	}
	
	public void setName(String name) {
		this.name = name;
		logger.info("Model parameter 'name' set to " + name);
	}
	
	public String getPathToCrux() {
		return pathToCrux;
	}
	
	public void setPathToCrux(String pathToCrux) {
		this.pathToCrux = pathToCrux;
	}
	
	private void setDefaults() {
		for (CruxComponents component: CruxComponents.values()) {
			componentsToRun[component.ordinal()] = false;
			showAdvancedParameters[component.ordinal()] = false;
			componentsRunStatus[component.ordinal()] = RunStatus.NOT_RUN;
		}
		allowMissedCleavages = allowMissedCleavagesDefault;
	    customEnzymeBeforeCleavage = customEnzymeBeforeCleavageDefault;
	    customEnzymeAfterCleavage = customEnzymeAfterCleavageDefault;
	    decoyPValues = decoyPValuesDefault;
	    digestType = digestTypeDefault;
	    enzyme = enzymeDefault;
	    featureFile = featureFileDefault;
	    massType = massTypeDefault;
	    maxLength = maxLengthDefault;
	    maxMass = maxMassDefault;
	    maxMods = maxModsDefaults;
	    minLength = minLengthDefault;
		minMass = minMassDefault;
		numDecoysPerTarget = numDecoysPerTargetDefault;
		numDecoyFiles = numDecoyFilesDefault;
	    pi0 = pi0Default;
	    printSearchProgress = printSearchProgressDefault;
		proteinSource = proteinSourceDefault;
		seed = seedDefault;
		spectraSource = spectraSourceDefault;
	    spectrumCharge = spectrumChargeDefault;
	    spectrumMaxMass = spectrumMaxMassDefault;
	    spectrumMinMass = spectrumMinMassDefault;
	    topMatch = topMatchDefault;
		logger.info("Model set to default values");
	}
	
	public void restoreDefaults(CruxComponents component) {
		componentsToRun[component.ordinal()] = false;
		showAdvancedParameters[component.ordinal()] = false;
		switch (component) {
		    case CREATE_INDEX:
        		allowMissedCleavages = allowMissedCleavagesDefault;
        	    customEnzymeBeforeCleavage = customEnzymeBeforeCleavageDefault;
	    		customEnzymeAfterCleavage = customEnzymeAfterCleavageDefault;
        		digestType = digestTypeDefault;
        		enzyme = enzymeDefault;
        		massType = massTypeDefault; 
	            maxLength = maxLengthDefault;
        	    maxMass = maxMassDefault;
	            minLength = minLengthDefault;
		        minMass = minMassDefault;
		        proteinSource = proteinSourceDefault;
		    	break;
		    case SEARCH_FOR_MATCHES:
        		allowMissedCleavages = allowMissedCleavagesDefault;
        	    customEnzymeBeforeCleavage = customEnzymeBeforeCleavageDefault;
	    		customEnzymeAfterCleavage = customEnzymeAfterCleavageDefault;
        	    decoyPValues = decoyPValuesDefault;
        		digestType = digestTypeDefault;
        		enzyme = enzymeDefault;
        		massType = massTypeDefault; 
	            maxLength = maxLengthDefault;
        	    maxMass = maxMassDefault;
	            maxMods = maxModsDefaults;
	            minLength = minLengthDefault;
		        minMass = minMassDefault;
        		numDecoysPerTarget = numDecoysPerTargetDefault;
				numDecoyFiles = numDecoyFilesDefault;
		        proteinSource = proteinSourceDefault;
		        spectraSource = spectraSourceDefault;
        	    spectrumCharge = spectrumChargeDefault;
        	    spectrumMaxMass = spectrumMaxMassDefault;
        	    spectrumMinMass = spectrumMinMassDefault;
		    	break;
		    case COMPUTE_Q_VALUES:
        	    pi0 = pi0Default;
		        proteinSource = proteinSourceDefault;
		    	break;
		    case PERCOLATOR:
        	    featureFile = featureFileDefault;
        	    pi0 = pi0Default;
		        proteinSource = proteinSourceDefault;
        	    topMatch = topMatchDefault;
		    	break;
		    case QRANKER:
        	    featureFile = featureFileDefault;
        	    pi0 = pi0Default;
		        proteinSource = proteinSourceDefault;
        	    topMatch = topMatchDefault;
		    	break;
		}
		logger.info("Model restored to default values");
	}
	
	public boolean getRunComponent(CruxComponents component) {
		return componentsToRun[component.ordinal()];
	}
	
	public void setRunComponent(CruxComponents component, boolean value) {
		componentsToRun[component.ordinal()] = value;
		if (value) {
		    ++numComponentsToRun;
		}
		else {
		    --numComponentsToRun;
		}
		logger.info("Run " + component.toString() + " set to " + value);
	}
	
	public boolean getShowAdvancedParameters(CruxComponents component) {
		return showAdvancedParameters[component.ordinal()];
	}
	
	public void setShowAdvancedParameters(CruxComponents component, boolean value) {
		showAdvancedParameters[component.ordinal()] = value;
		logger.info("Show advanced parameters for " + component.toString() + " set to " + value);
	}
	
	public IsotopicMassType getDefaultMassType() {
		return this.massTypeDefault;
	}
	
	public IsotopicMassType getMassType() {
		return this.massType;
	}

	public void setMassType(IsotopicMassType type) {
		this.massType = type;
		logger.info("Model parameter 'massType' set to " + type.toString());
	}
	
	public DigestType getDefaultDigestType() {
		return this.digestTypeDefault;
	}
	
	public DigestType getDigestType() {
		return this.digestType;
	}
	
	public void setDigestType(DigestType type) {
		this.digestType =type;
		logger.info("Model parameter 'digestType' set to " + type.toString());
	}
	
	public Enzyme getDefaultEnzyme() {
		return this.enzymeDefault;
	}
	
	public Enzyme getEnzyme() {
		return this.enzyme;
	}

	public void setEnzyme(Enzyme e) {
		this.enzyme = e;
		logger.info("Model parameter 'enzyme' set to " + e.toString());
	}
	
	public boolean getDefaultAllowMissedCleavages() {
		return this.allowMissedCleavagesDefault;
	}
	
	public boolean getAllowMissedCleavages() {
		return this.allowMissedCleavages;
	}
	
	public void setAllowMissedCleavages(boolean value) {
		this.allowMissedCleavages = value;
		logger.info("Model parameter 'allowMissedCleavages' set to " + value);
	}

	public String getDefaultProteinDatabase() {
		return proteinSourceDefault;
	}
	
	public String getProteinDatabase() {
		return proteinSource;
	}
	
	public void setProteinDatabase(String database) {
		this.proteinSource = database;
		logger.info("Model parameter 'proteinDatabase' set to " + database);
	}
	
	public String getDefaultSpectraFilename() {
		return spectraSourceDefault;
	}
	
	public String getSpectraFilename() {
		return spectraSource;
	}
	
	public void setSpectraFilename(String fileName) {
		this.spectraSource = fileName;
		logger.info("Model parameter 'spectraSource' set to " + fileName);
	}
	
	public String getCustomEnzymeAfterCleavage() {
		return customEnzymeAfterCleavage;
	}

	public void setCustomEnzymeAfterCleavage(String customEnzymeAfterCleavage) {
		this.customEnzymeAfterCleavage = customEnzymeAfterCleavage;
	}

	public String getCustomEnzymeBeforeCleavage() {
		return customEnzymeBeforeCleavage;
	}

	public void setCustomEnzymeBeforeCleavage(String customEnzymeBeforeCleavage) {
		this.customEnzymeBeforeCleavage = customEnzymeBeforeCleavage;
	}

	public boolean isDecoyPValues() {
		return decoyPValues;
	}

	public void setDecoyPValues(boolean decoyPValues) {
		this.decoyPValues = decoyPValues;
	}

	public Enzyme getEnzymeDefault() {
		return enzymeDefault;
	}

	public boolean isFeatureFile() {
		return featureFile;
	}

	public void setFeatureFile(boolean featureFile) {
		this.featureFile = featureFile;
	}

	public int getMaxLength() {
		return maxLength;
	}

	public void setMaxLength(int maxLength) {
		this.maxLength = maxLength;
	}

	public double getMaxMass() {
		return maxMass;
	}

	public void setMaxMass(double maxMass) {
		this.maxMass = maxMass;
	}

	public int getMaxMods() {
		return maxMods;
	}

	public void setMaxMods(int maxMods) {
		this.maxMods = maxMods;
	}

	public int getMinLength() {
		return minLength;
	}

	public void setMinLength(int minLength) {
		this.minLength = minLength;
	}

	public double getMinMass() {
		return minMass;
	}

	public void setMinMass(double minMass) {
		this.minMass = minMass;
	}

	public int getNumDecoyFiles() {
		return numDecoyFiles;
	}

	public void setNumDecoyFiles(int numDecoyFiles) {
		this.numDecoyFiles = numDecoyFiles;
	}

	public int getNumDecoysPerTarget() {
		return numDecoysPerTarget;
	}

	public void setNumDecoysPerTarget(int numDecoysPerTarget) {
		this.numDecoysPerTarget = numDecoysPerTarget;
	}

	public double getPi0() {
		return pi0;
	}

	public void setPi0(double pi0) {
		this.pi0 = pi0;
	}

	public int getPrintSearchProgress() {
		return printSearchProgress;
	}

	public void setPrintSearchProgress(int printSearchProgress) {
		this.printSearchProgress = printSearchProgress;
	}

	public double getSeed() {
		return seed;
	}

	public void setSeed(double seed) {
		this.seed = seed;
	}

	public AllowedSpectrumCharge getSpectrumCharge() {
		return spectrumCharge;
	}

	public void setSpectrumCharge(AllowedSpectrumCharge spectrumCharge) {
		this.spectrumCharge = spectrumCharge;
	}

	public double getSpectrumMaxMass() {
		return spectrumMaxMass;
	}

	public void setSpectrumMaxMass(double spectrumMaxMass) {
		this.spectrumMaxMass = spectrumMaxMass;
	}

	public double getSpectrumMinMass() {
		return spectrumMinMass;
	}

	public void setSpectrumMinMass(double spectrumMinMass) {
		this.spectrumMinMass = spectrumMinMass;
	}

	public int getTopMatch() {
		return topMatch;
	}

	public void setTopMatch(int topMatch) {
		this.topMatch = topMatch;
	}

	public String getCustomEnzymeAfterCleavageDefault() {
		return customEnzymeAfterCleavageDefault;
	}

	public String getCustomEnzymeBeforeCleavageDefault() {
		return customEnzymeBeforeCleavageDefault;
	}

	public DigestType getDigestTypeDefault() {
		return digestTypeDefault;
	}

	public IsotopicMassType getMassTypeDefault() {
		return massTypeDefault;
	}

	public int getMaxLengthDefault() {
		return maxLengthDefault;
	}

	public double getMaxMassDefault() {
		return maxMassDefault;
	}

	public int getMaxModsDefaults() {
		return maxModsDefaults;
	}

	public int getMinLengthDefault() {
		return minLengthDefault;
	}

	public double getMinMassDefault() {
		return minMassDefault;
	}

	public int getNumDecoyFilesDefault() {
		return numDecoyFilesDefault;
	}

	public int getNumDecoysPerTargetDefault() {
		return numDecoysPerTargetDefault;
	}

	public double getPi0Default() {
		return pi0Default;
	}

	public int getPrintSearchProgressDefault() {
		return printSearchProgressDefault;
	}

	public double getSeedDefault() {
		return seedDefault;
	}

	public AllowedSpectrumCharge getSpectrumChargeDefault() {
		return spectrumChargeDefault;
	}

	public double getSpectrumMaxMassDefault() {
		return spectrumMaxMassDefault;
	}

	public double getSpectrumMinMassDefault() {
		return spectrumMinMassDefault;
	}

	public int getTopMatchDefault() {
		return topMatchDefault;
	}
	
	public RunStatus getRunStatus(CruxComponents component) {
		return componentsRunStatus[component.ordinal()];
	}
	public void setRunStatus(CruxComponents component, RunStatus status) {
		componentsRunStatus[component.ordinal()] = status;
	}
}