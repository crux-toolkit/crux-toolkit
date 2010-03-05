package edu.washington.gs.noble.crux.gui;

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.util.logging.Logger;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JPanel;
import javax.swing.JTextField;

@SuppressWarnings("serial")
class SpectraPanel extends JPanel {
	
	JButton spectraButton = new JButton("Choose Spectra File");
	JTextField spectraFileName = new JTextField();
	private static Logger logger = 
		Logger.getLogger("edu.washington.gs.noble.crux.gui");
	
	class SpectraButtonListener implements ActionListener {
		public void actionPerformed(final ActionEvent event) {
			final JFileChooser chooser = new JFileChooser();
			chooser.setCurrentDirectory(new File("."));
			final int result = chooser.showOpenDialog(getParent());
			if (result == JFileChooser.APPROVE_OPTION){
				spectraFileName.setText(chooser.getSelectedFile().getPath());
				logger.info("User selected MS2 spectra file " + spectraFileName.getText() + ".");
			}
		}
	}
	
	public SpectraPanel(CruxAnalysisModel model) {
		setLayout(new BoxLayout(this, BoxLayout.LINE_AXIS));
		setBackground(Color.white);
		setAlignmentX(Component.LEFT_ALIGNMENT);
		spectraButton.setAlignmentX(Component.LEFT_ALIGNMENT);
		spectraButton.setAlignmentY(Component.CENTER_ALIGNMENT);
		spectraButton.addActionListener(new SpectraButtonListener());
		add(spectraButton);
		add(Box.createRigidArea(new Dimension(12,0)));
		spectraFileName.setAlignmentX(Component.LEFT_ALIGNMENT);
		spectraFileName.setAlignmentY(Component.CENTER_ALIGNMENT);
		add(spectraFileName);
		spectraFileName.setMaximumSize(new Dimension(200,25));
	}

}