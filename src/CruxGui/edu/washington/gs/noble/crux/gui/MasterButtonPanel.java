package edu.washington.gs.noble.crux.gui;

import java.awt.BorderLayout;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.util.logging.Logger;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import javax.swing.JPanel;

import edu.washington.gs.noble.crux.gui.CruxAnalysisModel.CruxComponents;

@SuppressWarnings("serial")
public class MasterButtonPanel extends JPanel {
	
	private static Logger logger = 
		Logger.getLogger("edu.washington.gs.noble.crux.gui");

	final CruxGui cruxGui;
	final JPanel buttonPanel = new JPanel();
	final JButton newButton = new JButton("New analysis");
	final JButton loadButton = new JButton("Load analysis");
	final JButton saveButton = new JButton("Save analysis");
	final JButton runButton = new JButton("Run analysis");
	final JButton setCruxPathButton = new JButton("Locate Crux");
	final JButton exitButton = new JButton("Exit");
	
	public MasterButtonPanel(final CruxGui cruxGui) {
		this.cruxGui = cruxGui;
		buttonPanel.setBorder(BorderFactory.createEtchedBorder());
		buttonPanel.setLayout(new FlowLayout());
		add(buttonPanel, BorderLayout.SOUTH);
		newButton.addActionListener(new NewAnalysisButtonListener()); 
		loadButton.addActionListener(new LoadAnalysisButtonListener()); 
		runButton.addActionListener(new RunButtonListener());
		saveButton.addActionListener(new SaveButtonListener()); 
		setCruxPathButton.addActionListener(new SetCruxPathListener());
		exitButton.addActionListener(new ExitButtonListener());
		buttonPanel.add(setCruxPathButton);
		buttonPanel.add(newButton);
		buttonPanel.add(saveButton);
		buttonPanel.add(loadButton);
		buttonPanel.add(runButton);
		buttonPanel.add(exitButton);
	}	
	
	class SetCruxPathListener implements ActionListener {
		
		public void actionPerformed(ActionEvent event) {
			cruxGui.promptForCruxPath();
		}
		
	}
	
	class NewAnalysisButtonListener implements ActionListener {

		public void actionPerformed(ActionEvent event) {
			CruxAnalysisModel model = cruxGui.getAnalysisModel();
    		boolean continueNewModel = true;
    		if (model != null && model.needsSaving()) {
    			continueNewModel = cruxGui.promptToSave();
    		}
    		if (continueNewModel) {
			    cruxGui.setAnalysisModel(new CruxAnalysisModel());
    		}
		}
		
	}
	
	class LoadAnalysisButtonListener implements ActionListener {

		public void actionPerformed(ActionEvent event) {
			CruxAnalysisModel model = cruxGui.getAnalysisModel();
    		boolean continueLoadModel = true;
    		if (model != null && model.needsSaving()) {
    			continueLoadModel = cruxGui.promptToSave();
    		}
    		if (continueLoadModel) {
    			String name = (String) JOptionPane.showInputDialog(
    					cruxGui.frame, 
    					"Choose the analysis to load:", 
    					"Load Analysis", 
    					JOptionPane.PLAIN_MESSAGE, 
    					null, 
    					cruxGui.modelSet.toArray(), 
    					""
    			);
    			if (name != null) {
        			model = cruxGui.getAnalysisModel().readModelFromBinaryFile(name);
        			if (model != null) {
        				cruxGui.setAnalysisModel(model);
        			}
    			}
    		}
		}
		
	}
	
	class SaveButtonListener implements ActionListener {

		public void actionPerformed(ActionEvent event) {
			boolean continueSave = true;
			CruxAnalysisModel model = cruxGui.getAnalysisModel();
        	// Model has to have name before we can save it
        	if (model.getName() == null) {
        		// Get name from user
        		continueSave = cruxGui.promptForModelName();
        	}
        	if (continueSave) {
        	    model.saveModelToBinaryFile();
        	}
		}
		
	}

	class RunButtonListener implements ActionListener {
		
		public void actionPerformed(ActionEvent event) {
			CruxAnalysisModel model = cruxGui.getAnalysisModel();
    		boolean continueRunAnalysis = true;
    		if (model != null && model.needsSaving()) {
    			continueRunAnalysis = cruxGui.promptToSave();
    		}
    		if (continueRunAnalysis) {
    			if (model.isValidAnalysis()) {
    				if (model.getPathToCrux() == null) {
    					continueRunAnalysis = cruxGui.promptForCruxPath();
    				}
    				if (continueRunAnalysis) {
    			        Process process = null;
						// Write out a parameter file to the analysis directory
						if (!model.saveModelToParameterFile()) {
							logger.info("Unable to write parameter file.");
							return;
						}
	
    			        for (CruxComponents component: CruxComponents.values()) {
	    			        process = model.runComponent(component);
	    			        if (process != null) {
		    			        try {
									process.waitFor();
								} catch (InterruptedException e) {
									// TODO Auto-generated catch block
									e.printStackTrace();
								}
								logger.info("Execution of " + component.toString() + " completed with exit status " + process.exitValue());
								if (process.exitValue() == 0) {
									model.setRunStatus(component, CruxAnalysisModel.RunStatus.COMPLETED);
								}
								else {
									model.setRunStatus(component, CruxAnalysisModel.RunStatus.FAILED);
									break;
								}
	    			        }
    			        }
    			        cruxGui.frame.updateFromModel(model);
    			        cruxGui.frame.repaint();
    				}
    			}
    		}
		}
		
	}
	
	class ExitButtonListener implements ActionListener {
		
		public void actionPerformed(final ActionEvent e) {
			cruxGui.shutdown();
		}
		
	}
}
