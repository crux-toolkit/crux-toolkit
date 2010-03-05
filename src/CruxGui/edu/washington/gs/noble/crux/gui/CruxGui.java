package edu.washington.gs.noble.crux.gui;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import java.util.logging.*;

public class CruxGui {

	private static Logger logger = Logger
			.getLogger("edu.washington.gs.noble.crux.gui");
	
	private static void shutdownGui() {
		logger.info("Shutting down");
		System.exit(0);
	}
	
	public static void main(String[] args) {

		logger.info("Staring CruxGui");
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				logger.info("Starting the Crux GUI Event Handler");
				CruxAnalysisModel model = new CruxAnalysisModel();
				JFrame frame = new JFrame();
				frame.setResizable(false);
				frame.setUndecorated(false);
				frame.addWindowListener(new WindowAdapter() {
					public void windowClosing(WindowEvent we) {
						shutdownGui();
					}
				});
				frame.setTitle("Crux");
				frame.getContentPane().setBackground(Color.white);
				frame.getContentPane().setLayout(new BorderLayout());
				CruxGuiPanel cruxGuiPanel = new  CruxGuiPanel(model);
				frame.add(cruxGuiPanel, BorderLayout.CENTER);
				JPanel controls = new JPanel();
				frame.add(controls, BorderLayout.SOUTH);
				controls.setSize(new Dimension(0, 40));
				//controls.setLayout(new BoxLayout(frame, BoxLayout.X_AXIS));
				JButton setName = new JButton("Set Analysis Name");
				setName.addActionListener(new NameButtonListener(model, frame)); 
				controls.add(setName);
				JCheckBox additionalParameters = new JCheckBox("Show Additional Parameters");
				controls.add(additionalParameters);
				controls.setBackground(additionalParameters.getBackground());
				JButton exitButton = new JButton("Exit");
				exitButton.addActionListener(new ActionListener() {
					public void actionPerformed(final ActionEvent e) {
						shutdownGui();
					}
				});
				controls.add(exitButton);
				frame.setSize(cruxGuiPanel.getWidth(), 4 * cruxGuiPanel.getHeight());
				frame.setLocationRelativeTo(null);
				frame.setVisible(true);
			}
		});
	}
	
}

class NameButtonListener implements ActionListener {

	CruxAnalysisModel model;
	JFrame frame;
	
	NameButtonListener(CruxAnalysisModel model, JFrame frame) {
		super();
		this.model = model;
		this.frame = frame;
	}
	
	public void actionPerformed(ActionEvent event) {
		String name = (String) JOptionPane.showInputDialog(
				frame, 
				"Choose a name for the analysis:", 
				"Analysis Name", 
				JOptionPane.PLAIN_MESSAGE, 
				null, 
				null, 
				""
		);
		model.setName(name);
		frame.setTitle("Crux - " + name);
	}	
}