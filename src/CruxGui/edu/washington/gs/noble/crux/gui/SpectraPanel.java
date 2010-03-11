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
	JTextField spectraFilename = new JTextField();
	CruxAnalysisModel model;
	
	private static Logger logger = 
		Logger.getLogger("edu.washington.gs.noble.crux.gui");
	
	class SpectraButtonListener implements ActionListener {
		public void actionPerformed(final ActionEvent event) {
			final JFileChooser chooser = new JFileChooser();
			String fileName = model.getSpectraFilemame();
			if (fileName != null && fileName.length() > 0) {
				chooser.setCurrentDirectory(new File(fileName));
			}
			else {
				chooser.setCurrentDirectory(new File("."));
			}
			final int result = chooser.showOpenDialog(getParent());
			if (result == JFileChooser.APPROVE_OPTION){
				spectraFilename.setText(chooser.getSelectedFile().getPath());
				logger.info("User selected MS2 spectra file " + spectraFilename.getText() + ".");
			}
		}
	}
	
	public SpectraPanel(CruxAnalysisModel model) {
		this.model = model;
		setLayout(new BoxLayout(this, BoxLayout.LINE_AXIS));
		setBackground(Color.white);
		setAlignmentX(Component.LEFT_ALIGNMENT);
		spectraButton.setAlignmentX(Component.LEFT_ALIGNMENT);
		spectraButton.setAlignmentY(Component.CENTER_ALIGNMENT);
		spectraButton.addActionListener(new SpectraButtonListener());
		add(spectraButton);
		add(Box.createRigidArea(new Dimension(12,0)));
		spectraFilename.setAlignmentX(Component.LEFT_ALIGNMENT);
		spectraFilename.setAlignmentY(Component.CENTER_ALIGNMENT);
		add(spectraFilename);
		spectraFilename.setMaximumSize(new Dimension(200,25));
	}

	public  String  getSpectraFilename() {
		return  spectraFilename.getText();
	}
	
	public void setSpectraFilename(String fileName) {
		spectraFilename.setText(fileName);
	}
	
	public void updateFromModel() {
		spectraFilename.setText(model.getSpectraFilemame());
	}

}