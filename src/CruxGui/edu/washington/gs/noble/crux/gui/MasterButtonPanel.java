package edu.washington.gs.noble.crux.gui;

import java.awt.BorderLayout;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JOptionPane;
import javax.swing.JPanel;

public class MasterButtonPanel extends JPanel {

	final JPanel buttonPanel = new JPanel();
	final JButton newButton = new JButton("New Analysis");
	final JButton openButton = new JButton("Open Analysis");
	final JButton setNameButton = new JButton("Set Analysis Name");
	final JButton runButton = new JButton("Run Analysis");
	final JButton exitButton = new JButton("Exit");
	CruxGui gui;
	
	public MasterButtonPanel(final CruxGui gui) {
		
		this.gui = gui;
		buttonPanel.setBorder(BorderFactory.createEtchedBorder());
		add(buttonPanel, BorderLayout.SOUTH);
		newButton.addActionListener(new OpenAnalysisButtonListener()); 
		openButton.addActionListener(new OpenAnalysisButtonListener()); 
		runButton.addActionListener(new RunButtonListener());
		setNameButton.addActionListener(new SetNameButtonListener()); 
		exitButton.addActionListener(new ActionListener() {
			public void actionPerformed(final ActionEvent e) {
				gui.shutdown();
			}
		});
		buttonPanel.setLayout(new FlowLayout());
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

	class SetNameButtonListener implements ActionListener {

		public void actionPerformed(ActionEvent event) {
			// Get list of existing names
			String name = (String) JOptionPane.showInputDialog(
					gui.frame, 
					"Choose a name for the analysis:", 
					"Analysis Name", 
					JOptionPane.PLAIN_MESSAGE, 
					null, 
					null, 
					""
			);
			if (name != null) {
				if (name == gui.model.getName()) {
					// Curren model alread has this name
					return;
				}
				if (gui.analysisSet.contains(name)) {
					JOptionPane.showMessageDialog(gui.frame, "An analysis named " + name + "already exists.");
					return;
				}
				// Create new directory using name
				boolean result = (new File(name)).mkdir();
				// Add name to analysis set and set name in model
				if (result) {
					gui.model.setName(name);
					gui.analysisSet.add(name);
					gui.frame.setTitle("Crux - " + name);
				}
				else {
					JOptionPane.showMessageDialog(gui.frame, "Unable to create directory for analysis named " + name);
				}
			}
		}
		
	}

	class RunButtonListener implements ActionListener {
		
		public void actionPerformed(ActionEvent event) {
			if (gui.model.isValidAnalysis()) {
			    gui.frame.cruxGuiPanel.setVisible(false);
			    gui.model.run();
			}
		}
		
	}
}
