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

public class MasterButtonPanel extends JPanel {

	final CruxGui cruxGui;
	final JPanel buttonPanel = new JPanel();
	final JButton newButton = new JButton("New Analysis");
	final JButton openButton = new JButton("Open Analysis");
	final JButton setNameButton = new JButton("Set Analysis Name");
	final JButton runButton = new JButton("Run Analysis");
	final JButton setCruxPathButton = new JButton("Locate Crux");
	final JButton exitButton = new JButton("Exit");
	
	public MasterButtonPanel(final CruxGui cruxGui) {
		this.cruxGui = cruxGui;
		buttonPanel.setBorder(BorderFactory.createEtchedBorder());
		add(buttonPanel, BorderLayout.SOUTH);
		newButton.addActionListener(new OpenAnalysisButtonListener()); 
		openButton.addActionListener(new OpenAnalysisButtonListener()); 
		runButton.addActionListener(new RunButtonListener());
		setNameButton.addActionListener(new SetNameButtonListener()); 
		setCruxPathButton.addActionListener(new SetCruxPathListener());
		exitButton.addActionListener(new ActionListener() {
			public void actionPerformed(final ActionEvent e) {
				cruxGui.shutdown();
			}
		});
		buttonPanel.setLayout(new FlowLayout());
		buttonPanel.add(setCruxPathButton);
		buttonPanel.add(setNameButton);
		buttonPanel.add(newButton);
		buttonPanel.add(openButton);
		buttonPanel.add(runButton);
		buttonPanel.add(exitButton);
	}	
	
	class OpenAnalysisButtonListener implements ActionListener {

		public void actionPerformed(ActionEvent event) {
		}
		
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
}
