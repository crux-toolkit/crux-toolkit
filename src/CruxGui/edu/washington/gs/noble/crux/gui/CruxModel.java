package edu.washington.gs.noble.crux.gui;

public class CruxModel {
	
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
	
	public enum DigestionType {
		FULL("full"), 
		PARTIAL("partial"); 
		String name;
		DigestionType(String name) {this.name = name;}
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
	}
	
	private boolean componentsToRun[]= new boolean[CruxComponents.getSize()];
	private String proteinSource;
	private String spectraSource;
	private IsotopicMassType massType;
	private DigestionType digestionType;
	private Enzyme enzyme;
	private boolean allowMissedCleavages;
	
	public CruxModel() {
		this.componentsToRun[CruxComponents.CREATE_INDEX.ordinal()] = false;
		this.componentsToRun[CruxComponents.SEARCH_FOR_MATCHES.ordinal()] = false;
		this.proteinSource = null;
		this.spectraSource = null;
		this.massType = IsotopicMassType.AVERAGE;
		this.digestionType = DigestionType.FULL;
		this.enzyme = Enzyme.TRYPSIN;
		this.allowMissedCleavages = false;
	}
	
	boolean getRunIndex() {
		return componentsToRun[CruxComponents.CREATE_INDEX.ordinal()];
	}

	void setRunIndex(boolean value) {
		componentsToRun[CruxComponents.CREATE_INDEX.ordinal()] = value;
	}
	
	boolean getRunSearch() {
		return componentsToRun[CruxComponents.SEARCH_FOR_MATCHES.ordinal()];
	}

	void setRunSearch(boolean value) {
		componentsToRun[CruxComponents.SEARCH_FOR_MATCHES.ordinal()] =value;
	}
	
	IsotopicMassType getMassType() {
		return this.massType;
	}

	void setMassType(IsotopicMassType type) {
		this.massType = type;
	}
	
	boolean getAllowMissedCleavages() {
		return this.allowMissedCleavages;
	}
	
	void setAllowMissedCleavages(boolean value) {
		this.allowMissedCleavages = value;
	}
}