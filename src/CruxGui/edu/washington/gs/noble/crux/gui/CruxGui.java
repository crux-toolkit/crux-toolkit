package edu.washington.gs.noble.crux.gui;

import java.awt.EventQueue;
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.logging.Logger;

/**
 * This is the main class for the CruxGui. The CruxGui is a GUI wrapper for the crux
 * set of tools for peptide-spectra matching. It allows the users to specify the most
 * common parameters for a crux analysis, sets up the parameter files and output directories,
 * and builds and executes the crux command lines.
 * 
 * The GUI interaction thread for the applications is started via the run() method for this class.
 *
 * @author Charles E. Grant
 *
 */

public class CruxGui implements Runnable {
	
	static final String ANALYIS_LIST_FILENAME = "crux-analyses.bin";
	public SortedSet<String> analysisSet;
	public CruxMainFrame frame;
	private CruxAnalysisModel model;

	private static Logger logger = Logger
			.getLogger("edu.washington.gs.noble.crux.gui");
	
	public static void main(String[] args) {
		EventQueue.invokeLater(new CruxGui());		
	}
		
	public void run() {
		logger.info("Starting the Crux GUI");
		boolean result = readModelSetFromFile();
		if (!result) {
			System.exit(-1);
		}
		model = new CruxAnalysisModel();
		frame = new CruxMainFrame(this);
		frame.setVisible(true);
	}
	
	public void shutdown() {
		logger.info("Shutting down the Crux GUI");
		System.exit(0);
	}
	
	public CruxAnalysisModel getAnalysisModel() {
		return model;
	}
	
	public void setAnalysisModel(CruxAnalysisModel model) {
		this.model = model;
	}
	
	@SuppressWarnings("unchecked")
	public boolean readModelSetFromFile() {
		logger.info("Reading list of existing analyses.");
		File nameListFile = new File(ANALYIS_LIST_FILENAME);
		if (nameListFile.exists()) {
			try {
				FileInputStream f = new FileInputStream(nameListFile);
		        ObjectInputStream o = new ObjectInputStream(f);
			    analysisSet = (TreeSet<String>) o.readObject();
			    o.close();
			    f.close();
			} catch(EOFException e1) {
				analysisSet = new TreeSet<String>();
			} catch (IOException e2) {
				logger.info("Unable to open file " 
						+ nameListFile.toString() 
						+ ": " + e2.toString());
				e2.printStackTrace();
				return false;
			} catch (ClassNotFoundException e3) {
				logger.info("Unable to open file " 
						+ nameListFile.toString() 
						+ ": " + e3.toString());
				e3.printStackTrace();
				return false;
		    }
		}
		else {
			analysisSet = new TreeSet<String>();
		}
		return true;
	}
	
	public boolean saveModelSetToFile() {
		
		logger.info("Reading list of existing analyses.");
		File nameListFile = new File(ANALYIS_LIST_FILENAME);
    	try {
    		FileOutputStream f = new FileOutputStream(nameListFile);
            ObjectOutputStream o = new ObjectOutputStream(f);
            o.writeObject(analysisSet);
            o.close();
            f.close();
    	} catch (IOException e1) {
    		logger.info("Unable to open file " 
    				+ nameListFile.toString() 
    				+ ": " + e1.toString());
    		e1.printStackTrace();
    		return false;
    	}
    	
    	return true;
	}
}


