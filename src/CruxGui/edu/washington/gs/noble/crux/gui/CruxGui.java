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

import javax.swing.JFileChooser;
import javax.swing.JOptionPane;

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
	public SortedSet<String> modelSet;
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
		frame.setFocusable(true);
		frame.setInitialFocus();
	}
	
	public boolean promptToSave() {
		
		boolean result = true;
		final String message = "You have unsaved changes, would you like to save them?";
		final String title = "Exiting with unsaved changes";
        int choice = JOptionPane.showOptionDialog(
        	frame, 
        	message, 
        	title, 
        	JOptionPane.YES_NO_CANCEL_OPTION, 
        	JOptionPane.QUESTION_MESSAGE, 
        	null, 
        	null, 
        	null
        );
        switch (choice) {
	        case JOptionPane.YES_OPTION:
	        	// Model has to have name before we can save it
	        	if (model.getName() == null) {
	        		// Get name from user
	        		result = promptForModelName();
	        	}
	        	if (result) {
	        	    model.saveModelToBinaryFile();
	        	}
	        	break;
	        case JOptionPane.CANCEL_OPTION:
	        	result = false;
	        	break;
	        case JOptionPane.NO_OPTION:
	        	break;
        }
        return result;
	}
	
	public boolean promptForCruxPath() {
		final JFileChooser fileChooser = new JFileChooser();
		fileChooser.setDialogTitle("Set path to crux executable");
	    int returnVal = fileChooser.showDialog(frame, "Set path to crux");
        if (returnVal == JFileChooser.APPROVE_OPTION) {
            File file = fileChooser.getSelectedFile();
            model.setPathToCrux(file.toString());
            logger.info("Set crux path to " + file.toString());
            return true;
        }
        else {
        	return false;
        }
	}
	
	public boolean promptForModelName() {
		
    	// Get name from user via dialog box
    	String name = (String) JOptionPane.showInputDialog(
    			frame, 
    			"Choose a name for the analysis:", 
    			"Analysis Name", 
    			JOptionPane.PLAIN_MESSAGE, 
    			null, 
    			null, 
    			""
    	);
    	
    	if (name == null) {
    		return false;
    	}
    	else {
    		
    		if (name == getAnalysisModel().getName()) {
    			// Current model already has this name, we don't need to do anything.
    			return true;
    		}
    		
    		if (modelSet.contains(name)) {
    			// An analysis with this name already exists, advise user
    			final String message = 
    				"An analysis named " + name + "already exists.\nPlease choose another name";
    			JOptionPane.showMessageDialog(frame, message);
    			return false;
    		}
    		
    		// Valid name
    		// Create directory for analysis using name
    		boolean result = (new File(name)).mkdir();
    		
    		if (result) {
    			//Get current model
    			CruxAnalysisModel model = getAnalysisModel();
    	
    			// Add name to analysis set and set name in model
    			model.setName(name);
    		    // File needs to be saved now
    			model.setNeedsSaving(true);
    			// Add model to list of known models.
    			modelSet.add(name);
    			// Save updated analysis set
    			saveModelSetToFile();
    			// Update frame title bar
    			frame.setTitle("Crux - " + name);
    			return true;
    		}
    		else {
    			String message = "Unable to create directory for analysis named " + name;
    			JOptionPane.showMessageDialog(frame, message);
    			return false;
    		}
    	}
	}

	public void shutdown() {
		boolean continueExit = true;
		if (model != null && model.needsSaving()) {
			continueExit = promptToSave();
		}
		if (continueExit) {
    		logger.info("Shutting down the Crux GUI");
    		System.exit(0);
		}
	}
	
	public CruxAnalysisModel getAnalysisModel() {
		return model;
	}
	
	public void setAnalysisModel(CruxAnalysisModel model) {
		this.model = model;
		frame.updateFromModel(model);
	}
	
	@SuppressWarnings("unchecked")
	public boolean readModelSetFromFile() {
		logger.info("Reading list of existing analyses.");
		File nameListFile = new File(ANALYIS_LIST_FILENAME);
		if (nameListFile.exists()) {
			try {
				FileInputStream f = new FileInputStream(nameListFile);
		        ObjectInputStream o = new ObjectInputStream(f);
			    modelSet = (TreeSet<String>) o.readObject();
			    o.close();
			    f.close();
			} catch(EOFException e1) {
				modelSet = new TreeSet<String>();
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
			modelSet = new TreeSet<String>();
		}
		return true;
	}
	
	public boolean saveModelSetToFile() {
		
		logger.info("Reading list of existing analyses.");
		File nameListFile = new File(ANALYIS_LIST_FILENAME);
    	try {
    		FileOutputStream f = new FileOutputStream(nameListFile);
            ObjectOutputStream o = new ObjectOutputStream(f);
            o.writeObject(modelSet);
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


