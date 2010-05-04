package edu.washington.gs.noble.crux.gui;

import java.awt.BorderLayout;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import javax.swing.JPanel;

@SuppressWarnings("serial")
public class MasterButtonPanel extends JPanel {

	final CruxGui cruxGui;
	final JPanel buttonPanel = new JPanel();
	final JButton newButton = new JButton("New Analysis");
	final JButton loadButton = new JButton("Load Analysis");
	final JButton setNameButton = new JButton("Set Analysis Name");
	final JButton runButton = new JButton("Run Analysis");
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
		setNameButton.addActionListener(new SetNameButtonListener()); 
		setCruxPathButton.addActionListener(new SetCruxPathListener());
		exitButton.addActionListener(new ExitButtonListener());
		buttonPanel.add(setCruxPathButton);
		buttonPanel.add(newButton);
		buttonPanel.add(setNameButton);
		buttonPanel.add(loadButton);
		buttonPanel.add(runButton);
		buttonPanel.add(exitButton);
	}	
	
	class SetCruxPathListener implements ActionListener {
		
		public void actionPerformed(ActionEvent event) {
			final JFileChooser fileChooser = new JFileChooser();
			    int returnVal = fileChooser.showDialog(getParent(), "Set path to crux");
		        if (returnVal == JFileChooser.APPROVE_OPTION) {
		            File file = fileChooser.getSelectedFile();
		             cruxGui.getAnalysisModel().setPathToCrux(file.toString());
		        }
		}
		
	}
	
	class NewAnalysisButtonListener implements ActionListener {

		public void actionPerformed(ActionEvent event) {
			cruxGui.setAnalysisModel(new CruxAnalysisModel());
		}
		
	}
	
	class LoadAnalysisButtonListener implements ActionListener {

		public void actionPerformed(ActionEvent event) {
			String name = (String) JOptionPane.showInputDialog(
					cruxGui.frame, 
					"Choose the analysis to load:", 
					"Load Analysis", 
					JOptionPane.PLAIN_MESSAGE, 
					null, 
					cruxGui.analysisSet.toArray(), 
					""
			);
			CruxAnalysisModel model = cruxGui.getAnalysisModel().readModelFromBinaryFile(name);
			if (model != null) {
				cruxGui.setAnalysisModel(model);
			}
		}
		
	}
	
	class SetNameButtonListener implements ActionListener {

		public void actionPerformed(ActionEvent event) {
			// Get list of existing names
			String name = (String) JOptionPane.showInputDialog(
					cruxGui.frame, 
					"Choose a name for the analysis:", 
					"Analysis Name", 
					JOptionPane.PLAIN_MESSAGE, 
					null, 
					null, 
					""
			);
			if (name != null) {
				if (name == cruxGui.getAnalysisModel().getName()) {
					// Current model already has this name
					return;
				}
				if (cruxGui.analysisSet.contains(name)) {
					JOptionPane.showMessageDialog(cruxGui.frame, "An analysis named " + name + "already exists.");
					return;
				}
				// Create new directory using name
				boolean result = (new File(name)).mkdir();
				// Add name to analysis set and set name in model
				if (result) {
					cruxGui.getAnalysisModel().setName(name);
					cruxGui.analysisSet.add(name);
					cruxGui.frame.setTitle("Crux - " + name);
				}
				else {
					JOptionPane.showMessageDialog(cruxGui.frame, "Unable to create directory for analysis named " + name);
				}
			}
		}
		
	}

	class RunButtonListener implements ActionListener {
		
		public void actionPerformed(ActionEvent event) {
			if (cruxGui.getAnalysisModel().isValidAnalysis()) {
			    cruxGui.getAnalysisModel().run();
			}
		}
		
	}
	
	class ExitButtonListener implements ActionListener {
		
			public void actionPerformed(final ActionEvent e) {
				cruxGui.shutdown();
			}
		
	}
}
