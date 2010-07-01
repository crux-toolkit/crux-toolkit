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
	final JButton setupButton = new JButton("Set up analysis");
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
		setupButton.addActionListener(new SetupButtonListener());
		setupButton.setVisible(false);
		saveButton.addActionListener(new SaveButtonListener()); 
		setCruxPathButton.addActionListener(new SetCruxPathListener());
		exitButton.addActionListener(new ExitButtonListener());
		buttonPanel.add(setCruxPathButton);
		buttonPanel.add(newButton);
		buttonPanel.add(saveButton);
		buttonPanel.add(loadButton);
		buttonPanel.add(runButton);
		buttonPanel.add(setupButton);
		buttonPanel.add(exitButton);
	}
    public void setInitialFocus(){
	setCruxPathButton.requestFocusInWindow();
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
	
	public void showRunAnalysis() {
		cruxGui.frame.cruxRunPanel.setVisible(true);
		cruxGui.frame.add(cruxGui.frame.cruxRunPanel, BorderLayout.CENTER);
		cruxGui.frame.cruxSetupPanel.setVisible(false);
		setCruxPathButton.setVisible(false);
		newButton.setVisible(false);
		saveButton.setVisible(false);
		loadButton.setVisible(false);
		runButton.setVisible(false);
		setupButton.setVisible(true);
		repaint();
	}
	
	public void showSetupAnalysis() {
		cruxGui.frame.cruxSetupPanel.setVisible(true);
		cruxGui.frame.add(cruxGui.frame.cruxSetupPanel, BorderLayout.CENTER);
		cruxGui.frame.cruxRunPanel.setVisible(false);
		setCruxPathButton.setVisible(true);
		newButton.setVisible(true);
		saveButton.setVisible(true);
		loadButton.setVisible(true);
		runButton.setVisible(true);
		setupButton.setVisible(false);
		repaint();
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
			if (model.isValidAnalysis()) {
				showRunAnalysis();
				logger.info("Starting Analysis.");
				cruxGui.frame.cruxRunPanel.runAnalysis();
				logger.info("Analysis completed.");
			}
		}
		
	}
	
	class SetupButtonListener implements ActionListener {
		
		public void actionPerformed(ActionEvent event) {
			logger.info("Switching to analysis setup");
			showSetupAnalysis();
		}
		
	}
	
	class ExitButtonListener implements ActionListener {
		
		public void actionPerformed(final ActionEvent e) {
			cruxGui.shutdown();
		}
		
	}
}
