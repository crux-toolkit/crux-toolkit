package edu.washington.gs.noble.crux.gui;

import java.io.Serializable;
import java.util.logging.Logger;

public class CruxAnalysisModel extends Object implements Serializable{
	
	private static final long serialVersionUID = 3352659336939985777L;
	private static Logger logger = 
		Logger.getLogger("edu.washington.gs.noble.crux.gui");


	public enum CruxComponents {
		CREATE_INDEX("create-index"), 
		SEARCH_FOR_MATCHES("search-for-matches");
		String name;
		CruxComponents(String name) {this.name = name;}
		public static int getSize() { return 2; };
	};
	
	public enum IsotopicMassType { 
		AVERAGE("average"), 
		MONOISOTOPIC("monoisotopic");
		String name;
		IsotopicMassType(String name) {this.name = name;}
		public static int getSize() { return 2; };
	};
	
	public enum DigestType {
		FULL("full"), 
		PARTIAL("partial"); 
		String name;
		DigestType(String name) {this.name = name;}
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
	
	void toCruxParameterFile() {
		
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
		logger.info("Model parameter 'IsotopicMassType' set to " + type.toString());
	}
	
	DigestType getDigestType() {
		return this.digestType;
	}
	
	void setDigestType(DigestType type) {
		this.digestType =type;
		logger.info("Model parameter 'DigestType' set to " + type.toString());
	}
	
	Enzyme getEnzyme() {
		return this.enzyme;
	}

	void setEnzyme(Enzyme e) {
		this.enzyme = e;
		logger.info("Model parameter 'Enzyme' set to " + e.toString());
	}
	
	boolean getAllowMissedCleavages() {
		return this.allowMissedCleavages;
	}
	
	void setAllowMissedCleavages(boolean value) {
		this.allowMissedCleavages = value;
		logger.info("Model parameter 'AllowMissedCleavages' set to " + value);
	}

}