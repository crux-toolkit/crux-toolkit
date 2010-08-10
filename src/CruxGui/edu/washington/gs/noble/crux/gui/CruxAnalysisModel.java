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
import java.util.HashSet;

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
	
	/** The list of amino acids */
	public enum AminoAcid{
		LYSINE("K"),
		ASPARAGINE("N"),
		THREONINE("T"),
		ARGININE("R"),
		METHIONINE("M"),
		ISOLEUCINE("I"),
		GLUTAMINE("Q"),
		HISTIDINE("H"),
		PROLINE("P"),
		GLUTAMIC("E"),
		ASPARTIC("D"),
		ALANINE("A"),
		GLYCINE("G"),
		VALINE("V"),
		TYROSINE("Y"),
		SERINE("S"),
		TRYPTOPHAN("W"),
		CYSTEINE("C"),
		LEUCINE("L"),
		PHENYLALANINE("F");
		String name;
		AminoAcid(String name) {this.name = name;};
		public static int getSize() {return 20;};
		public String toString() {return this.name;}
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
	
	/** The list of allowed spectrum charges */
	public enum DecoyLocation {
		TARGET_FILE("target-file"), 
		ONE_DECOY_FILE("one-decoy-file"),
		SEPARATE_DECOY_FILES("separate-decoy-files") ; 
		String name;
		DecoyLocation(String name) {this.name = name;}
		public String toString() { return this.name; }
		public static int getSize() { return 3; };
	};
    



    /** The list of allowed verbosity levels */
    public enum Verbosity{
    	FATAL_ERRORS("fatal errors"),
	    NON_FATAL_ERRORS("non-fatal errors"),
	    WARNINGS("warnings"),
	    INFORMATION_ON_THE_PROGRESS_OF_EXECUTION("info on the progress of execution"),
	    MORE_PROGRESS_INFORMATION("more progress information"),
	    DEBUG_INFO("debug info"),
	    DETAILED_DEBUG_INFO("detailed debug info");
    	String name;
	Verbosity(String name) {this.name = name;};
	public String toString() {return this.name; }
	public static int getSize() {return 6;}
	public String getLevel() {
	    if (this.name.equals("fatal errors")) return "0";
	    if (this.name.equals("non-fatal errors")) return "10";
	    if (this.name.equals("warnings")) return "20";
	    if (this.name.equals("info on the progress of execution")) return "30";
	    if (this.name.equals("more progress information")) return "40";
	    if (this.name.equals("debug info")) return "50";
	    if (this.name.equals("detailed debug info")) return "60";
	    return "30";
	}
	    
    };
	
	private boolean needsSaving = false;
	private String name;
	private String pathToCrux;
	private final boolean componentsToRun[]= new boolean[CruxComponents.getSize()];
	private final RunStatus componentsRunStatus[]= new RunStatus[CruxComponents.getSize()];
	private final boolean showAdvancedParameters[]= new boolean[CruxComponents.getSize()];

	
	// Crux parameters
	private HashSet<AminoAcid> selectedCustomEnzymeBefore;
	private boolean preventCustomEnzymeBefore;
	private HashSet<AminoAcid> selectedCustomEnzymeAfter;
	private boolean preventCustomEnzymeAfter;
	private boolean enableCustomEnzyme;
	private boolean allowMissedCleavages;
	private String outputDir;
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
	private DecoyLocation decoyLocation;
	private double pi0;
	private int printSearchProgress;
	private String proteinSource;
	private double seed;
	private String spectraSource;
	private AllowedSpectrumCharge spectrumCharge;
	private double spectrumMaxMass;
	private double spectrumMinMass;
	private int topMatch;
	private Verbosity verbosity;
        
	
	// Crux parameters defaults
	final private HashSet<CruxAnalysisModel.AminoAcid> selectedCustomEnzymeBeforeDefault = 
		new HashSet<CruxAnalysisModel.AminoAcid>();
	final private boolean preventCustomEnzymeBeforeDefault = false;
	final private HashSet<CruxAnalysisModel.AminoAcid> selectedCustomEnzymeAfterDefault =
		new HashSet<CruxAnalysisModel.AminoAcid>();
	final private boolean preventCustomEnzymeAfterDefault = false;
	final private boolean enableCustomEnzymeDefault = false;
	final private boolean allowMissedCleavagesDefault = false;
	final private String outputDirDefault = null;
	final private DecoyLocation decoyLocationDefault = DecoyLocation.SEPARATE_DECOY_FILES;
	final private boolean decoyPValuesDefault = false;
	final private DigestType digestTypeDefault = DigestType.FULL;
	final private Enzyme enzymeDefault = Enzyme.TRYPSIN;
	final private boolean featureFileDefault = false;
	final private IsotopicMassType massTypeDefault = IsotopicMassType.AVERAGE;
	final private int maxLengthDefault = 50;
	final private double maxMassDefault = 7200.0;
	final private int maxModsDefaults = 255;
	final private int minLengthDefault = 6;
	final private double minMassDefault = 200;
	final private int numDecoysPerTargetDefault = 2;
	final private int printSearchProgressDefault = 10;
    final private String proteinSourceDefault = null;
	final private String spectraSourceDefault = null;
	final private AllowedSpectrumCharge spectrumChargeDefault = AllowedSpectrumCharge.ALL;
	final private double seedDefault = -1.0;
	final private double spectrumMaxMassDefault = 1000000000.000000;
	final private double spectrumMinMassDefault = 0.0;
	final private int topMatchDefault = 5;
	final private Verbosity verbosityDefault=Verbosity.WARNINGS;
 
	
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

	public boolean saveModelToParameterFile() {
		boolean result = false;
		String fileName = name + "/" + "analysis.params";
		try {
			PrintWriter out = new PrintWriter(new File(fileName));
			out.println("#Parameter file for analysis " + name);
			if (allowMissedCleavages) {
				out.println("missed-cleavages=T");
			}
			else {
				out.println("missed-cleavages=F");
			}
			out.println("digestion=" + digestType);
			out.println("enzyme=" + enzyme);
			if (featureFile) {
				out.println("feature-file=T");
			}
			else {
				out.println("feature-file=F");
			}
			out.println("isotopic-mass=" + massType);
			out.println("max-length=" + maxLength);
			out.println("max-mass=" + maxMass);
			out.println("max-mods=" + maxMods);
			out.println("min-length=" + minLength);
			out.println("min-mass=" + minMass);
			out.println("num-decoys-per-target=" + numDecoysPerTarget);
			out.println("decoy-location=" + decoyLocation.toString());
			out.println("print-search-progress=" + printSearchProgress);
			out.println("spectrum-charge=" + spectrumCharge.toString());
			out.println("spectrum-max-mass=" + spectrumMaxMass);
			out.println("spectrum-min-mass=" + spectrumMinMass);
			out.println("top-match=" + topMatch);
			if (enableCustomEnzyme){
				out.println("custom-enzyme=" + generateCustomEnzyme());
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
		boolean indexExists = (new File(name +"/index")).exists();
		boolean result = true;
		
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
			if (spectraSource == null) {
				JOptionPane.showMessageDialog(null, "A spectra source must be specified to search for matches.");
				result = false;
			}

			if (!componentsToRun[CruxComponents.CREATE_INDEX.ordinal()] && !indexExists){
				JOptionPane.showMessageDialog(null, "The create index step must be run before search for matches.");
				result = false;
			}
		}
		if (componentsToRun[CruxComponents.COMPUTE_Q_VALUES.ordinal()]){
		    if (!componentsToRun[CruxComponents.CREATE_INDEX.ordinal()] && !indexExists){
				JOptionPane.showMessageDialog(null, "The create index step must be run before running compute q values.");
				result = false;
			}
		}
		if (componentsToRun[CruxComponents.PERCOLATOR.ordinal()]){
		    if (!componentsToRun[CruxComponents.CREATE_INDEX.ordinal()] && !indexExists){
				JOptionPane.showMessageDialog(null, "The create index step must be run before running percolator.");
				result = false;
			}
		}
		if (componentsToRun[CruxComponents.QRANKER.ordinal()]){
		    if (!componentsToRun[CruxComponents.CREATE_INDEX.ordinal()] && !indexExists){
				JOptionPane.showMessageDialog(null, "The create index step must be run before running q-ranker.");
				result = false;
			}
		}
			
		logger.info("Validation of analysis: " + result + ".");

		return result;
	}
	
	//  Run the component using the component directory and parameter file
	Process runComponent(CruxComponents component) {
		
		Process process = null;
		if (componentsToRun[component.ordinal()]  && componentsRunStatus[component.ordinal()] != RunStatus.COMPLETED) {
			String subCommand = component.toString();
			String command = null;
			String proteinSource = this.proteinSource;
			if (component != CruxComponents.CREATE_INDEX && componentsToRun[CruxComponents.CREATE_INDEX.ordinal()]) {
				// If we have an index and this is not the create index command use the
				// index for the protein source.
				if (componentsRunStatus[CruxComponents.CREATE_INDEX.ordinal()] == RunStatus.COMPLETED) {
					proteinSource = name + "/index";
				}
				else {
					logger.info("Index creation failed. Unable to execute command: " + component.toString());
					JOptionPane.showMessageDialog(null, "Index creation failed. Unable to execute command: " + component.toString());
					return process;
				}
			}
			switch(component) {
			         case CREATE_INDEX:
				     proteinSource = this.proteinSource;
				     command = pathToCrux 
						+ " " + subCommand 
					    + " --overwrite T " 
					    + " --verbosity "+ verbosity.getLevel()
						+ " --parameter-file " + name + "/analysis.params " 
						+ proteinSource 
						+ " " + name + "/index";
					break;
				case SEARCH_FOR_MATCHES:
					command = pathToCrux 
						+ " " + subCommand 
					 	+ " --overwrite T "
					    + " --verbosity " + verbosity.getLevel()
						+ " --output-dir " + name + "/crux-output"
					 	+ " --parameter-file " + name + "/analysis.params" 
					 	+ " " +spectraSource 
					 	+ " " +proteinSource;
					break;
				case COMPUTE_Q_VALUES:
					command = pathToCrux
						+ " " +subCommand
						+ " --overwrite T "
						+ " --verbosity " + verbosity.getLevel()
						+ " --output-dir " + name + "/crux-output"
						+ " --parameter-file " + name + "/analysis.params "
						+ " " +proteinSource
						+ " " +name + "/crux-output";
					break;
				case PERCOLATOR:
					command = pathToCrux
						+ " "+subCommand
						+ " --overwrite T "
						+ " --verbosity " + verbosity.getLevel()
						+ " --output-dir " + name + "/crnext step: gather relavant statistics on peak area calculationux-output"
						+ " --parameter-file " + name + "/analysis.params "
						+ " " +proteinSource
						+ " " +name+ "/crux-output";
					break;
				case QRANKER:
					command = pathToCrux 
						+ " " + subCommand 
					 	+ " --overwrite T "
					 	+ " --verbosity "+ verbosity.getLevel()
						+ " --output-dir " + name + "/crux-output"
					 	+ " --parameter-file " + name + "/analysis.params " 
					 	+ " " +proteinSource
					 	+ " " +name+ "/crux-output";
					break;
			}
			try {
				logger.info("Attempt to execute command: " + command);
				process = Runtime.getRuntime().exec(command);
			} catch (Exception e) {
				process = null;
				logger.info("Unable to execute command: " + e.toString());
				JOptionPane.showMessageDialog(null, "Unable to execute command: " + e.toString());
			}
		}
		return process;
		
	}
	
	private String generateCustomEnzyme(){

		String beforeOpenBracket = "[";
		String afterOpenBracket = "[";
		String beforeCloseBracket = "]";
		String afterCloseBracket = "]";
		if (preventCustomEnzymeBefore){
			beforeOpenBracket = "{";
			beforeCloseBracket = "}";
		} 
		if (preventCustomEnzymeAfter){
			afterOpenBracket = "[";
			afterCloseBracket = "]";
		}
		String beforeAmino = new String();
		for (AminoAcid aminoAcid : selectedCustomEnzymeBefore){
			beforeAmino += aminoAcid.toString();
		}
		String afterAmino = new String();
		for (AminoAcid aminoAcid : selectedCustomEnzymeAfter){
			afterAmino += aminoAcid.toString();
		}
		return beforeOpenBracket+beforeAmino+beforeCloseBracket+"|"+
			afterOpenBracket+afterAmino+afterCloseBracket;
		
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
		selectedCustomEnzymeBefore = selectedCustomEnzymeBeforeDefault;
		preventCustomEnzymeBefore = preventCustomEnzymeBeforeDefault;
		selectedCustomEnzymeAfter = selectedCustomEnzymeAfterDefault;
		preventCustomEnzymeAfter = preventCustomEnzymeAfterDefault;
		enableCustomEnzyme = enableCustomEnzymeDefault;
		allowMissedCleavages = allowMissedCleavagesDefault;
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
		decoyLocation = decoyLocationDefault;
	    printSearchProgress = printSearchProgressDefault;
	    proteinSource = proteinSourceDefault;
		seed = seedDefault;
		spectraSource = spectraSourceDefault;
	    spectrumCharge = spectrumChargeDefault;
	    spectrumMaxMass = spectrumMaxMassDefault;
	    spectrumMinMass = spectrumMinMassDefault;
	    topMatch = topMatchDefault;
	    verbosity = verbosityDefault;
		logger.info("Model set to default values");
	}
	
	public void restoreDefaults(CruxComponents component) {
		componentsToRun[component.ordinal()] = false;
		showAdvancedParameters[component.ordinal()] = false;
		verbosity = verbosityDefault;
		switch (component) {
		    case CREATE_INDEX:
		    	selectedCustomEnzymeBefore = selectedCustomEnzymeBeforeDefault;
				preventCustomEnzymeBefore = preventCustomEnzymeBeforeDefault;
				selectedCustomEnzymeAfter = selectedCustomEnzymeAfterDefault;
				preventCustomEnzymeAfter = preventCustomEnzymeAfterDefault;
				enableCustomEnzyme = enableCustomEnzymeDefault;
		    	allowMissedCleavages = allowMissedCleavagesDefault;
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
		    	decoyPValues = decoyPValuesDefault;
		    	maxMods = maxModsDefaults;
		    	numDecoysPerTarget = numDecoysPerTargetDefault;
		    	decoyLocation = decoyLocationDefault;
		    	spectraSource = spectraSourceDefault;
		    	spectrumCharge = spectrumChargeDefault;
		    	spectrumMaxMass = spectrumMaxMassDefault;
		    	spectrumMinMass = spectrumMinMassDefault;
		    	break;
		    case COMPUTE_Q_VALUES:
		    	break;
		    case PERCOLATOR:
        	    featureFile = featureFileDefault;
        	    topMatch = topMatchDefault;
		    	break;
		    case QRANKER:
        	    featureFile = featureFileDefault;
        	    topMatch = topMatchDefault;
		    	break;
		}
		logger.info("Model restored to default values");
	}
	
	
	
	
	
	
	
/****************************************************/
/*              Getters and Setters                 */
/****************************************************/
	
	public void setCustomEnzymeCleavage(HashSet<AminoAcid> selectedBefore, 
		HashSet<AminoAcid> selectedAfter, boolean preventBefore, boolean preventAfter, boolean enable){
		selectedCustomEnzymeBefore = selectedBefore;
		preventCustomEnzymeBefore = preventBefore;
		selectedCustomEnzymeAfter = selectedAfter;
		preventCustomEnzymeAfter = preventAfter;
		enableCustomEnzyme = enable;
	}
	
	public HashSet<AminoAcid>  getSelectedCustomEnzymeBefore(){
		return selectedCustomEnzymeBefore;
	}
	
	public HashSet<AminoAcid>  getSelectedCustomEnzymeAfter(){
		return selectedCustomEnzymeAfter;
	}
	
	public boolean getPreventCustomEnzymeBefore(){
		return preventCustomEnzymeBefore;
	}
	
	public boolean getPreventCustomEnzymeAfter(){
		return preventCustomEnzymeAfter;
	}
	
	public boolean getEnableCustomEnzyme(){
		return enableCustomEnzyme;
	}
	
	public HashSet<AminoAcid>  getDefaultSelectedCustomEnzymeBefore(){
		return selectedCustomEnzymeBeforeDefault;
	}
	
	public HashSet<AminoAcid>  getDefaultSelectedCustomEnzymeAfter(){
		return selectedCustomEnzymeAfterDefault;
	}
	
	public boolean getDefaultPreventCustomEnzymeBefore(){
		return preventCustomEnzymeBeforeDefault;
	}
	
	public boolean getDefaultPreventCustomEnzymeAfter(){
		return preventCustomEnzymeAfterDefault;
	}
	
	public boolean getDefaultEnableCustomEnzyme(){
		return enableCustomEnzymeDefault;
	}
	
	public boolean getRunComponent(CruxComponents component) {
		return componentsToRun[component.ordinal()];
	}
	
	public void setRunComponent(CruxComponents component, boolean value) {
		componentsToRun[component.ordinal()] = value;
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
	
	public boolean isDecoyPValues() {
		return decoyPValues;
	}

	public void setDecoyPValues(boolean decoyPValues) {
		this.decoyPValues = decoyPValues;
	}

	public Enzyme getEnzymeDefault() {
		return enzymeDefault;
	}

	public boolean useFeatureFile() {
		return featureFile;
	}

	public boolean getFeatureFileDefault() {
		return featureFileDefault;
	}

	public void setUseFeatureFile(boolean featureFile) {
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

	public DecoyLocation getDecoyLocation() {
		return decoyLocation;
	}

	public void setDecoyLocation(DecoyLocation location) {
		this.decoyLocation = location;
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

	public void setVerbosity(Verbosity verbosity){
		this.verbosity = verbosity;
	}

	public Verbosity getVerbosity(){
		return verbosity;
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


        public String getOutputDir(){
	        return outputDir;
	}

        public void setOutputDir(String outputDir){
	        this.outputDir = outputDir;
	}


	public DigestType getDigestTypeDefault() {
		return digestTypeDefault;
	}

	public IsotopicMassType getMassTypeDefault() {
		return massTypeDefault;
	}

	public String getOutputDirDefault(){
		return outputDirDefault;
	}
    
	public int getMaxLengthDefault() {
		return maxLengthDefault;
	}

	public double getMaxMassDefault() {
		return maxMassDefault;
	}

	public int getMaxModsDefault() {
		return maxModsDefaults;
	}

	public int getMinLengthDefault() {
		return minLengthDefault;
	}

	public double getMinMassDefault() {
		return minMassDefault;
	}

	public DecoyLocation getDecoyLocationDefault() {
		return decoyLocationDefault;
	}

	public int getNumDecoysPerTargetDefault() {
		return numDecoysPerTargetDefault;
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
	
	public Verbosity getVerbosityDefault(){
    	return verbosityDefault;
    }
    
	public RunStatus getRunStatus(CruxComponents component) {
		return componentsRunStatus[component.ordinal()];
	}
	public void setRunStatus(CruxComponents component, RunStatus status) {
		componentsRunStatus[component.ordinal()] = status;
	}


}