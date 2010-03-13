package edu.washington.gs.noble.crux.gui;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.util.logging.Logger;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;

@SuppressWarnings("serial")
public class CruxMainFrame extends JFrame {

	private static Logger logger = Logger
		.getLogger("edu.washington.gs.noble.crux.gui");

	CruxAnalysisModel model;
	CruxViewPanel cruxGuiPanel;
	
	public  CruxMainFrame(final CruxGui gui) {
		setResizable(false);
		setUndecorated(false);
		setTitle("Crux");
		setBackground(Color.white);
		setLayout(new BorderLayout());	
		addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent we) {
				gui.shutdown();
			}
		});
		model = gui.getAnalysisModel();
		cruxGuiPanel = new CruxViewPanel(model);
		add(cruxGuiPanel, BorderLayout.CENTER);
		JPanel controls = new JPanel();
		controls.setBorder(BorderFactory.createEtchedBorder());
		add(controls, BorderLayout.SOUTH);
		controls.setLayout(new BorderLayout());
		JButton setName = new JButton("Set Analysis Name");
		setName.addActionListener(new NameButtonListener()); 
		controls.add(setName, BorderLayout.WEST);
		JButton runAnalysis = new JButton("Run Analysis");
		runAnalysis.addActionListener(new RunButtonListener());
		controls.add(runAnalysis, BorderLayout.CENTER);
		JButton exitButton = new JButton("Exit");
		exitButton.addActionListener(new ActionListener() {
			public void actionPerformed(final ActionEvent e) {
				gui.shutdown();
			}
		});
		controls.add(exitButton, BorderLayout.EAST);
		setSize(cruxGuiPanel.getWidth(), 4 * cruxGuiPanel.getHeight());
		setLocationRelativeTo(null);

	}
	
	class NameButtonListener implements ActionListener {

		public void actionPerformed(ActionEvent event) {
			// Get list of existing names
			String name = (String) JOptionPane.showInputDialog(
					CruxMainFrame.this, 
					"Choose a name for the analysis:", 
					"Analysis Name", 
					JOptionPane.PLAIN_MESSAGE, 
					null, 
					null, 
					""
			);
			if (name != null) {
				// If name is in list of existing names
				// load existing analysis
				// otherwise create new directory using name
				// and set name in model
				boolean result = false;
				try {
					result = new File(name).mkdir();
				}
				catch (SecurityException e) {
					// Show message box and return
					JOptionPane.showMessageDialog(CruxMainFrame.this, "Unable to create directory for analysis name " + name + ": " + e.toString());
					return;
				}
				if (result) {
					model.setName(name);
					CruxMainFrame.this.setTitle("Crux - " + name);
				}
				else {
					JOptionPane.showMessageDialog(CruxMainFrame.this, "Unable to create directory for analysis name " + name);
				}
			}
		}
		
	}

	class RunButtonListener implements ActionListener {
		
		public void actionPerformed(ActionEvent event) {
			if (model.isValidAnalysis()) {
			    cruxGuiPanel.setVisible(false);
			    model.run();
			}
		}
		
	}

}
;