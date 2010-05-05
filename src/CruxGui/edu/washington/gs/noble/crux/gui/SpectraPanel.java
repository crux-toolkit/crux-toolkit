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
class SpectraPanel extends CruxParameterControl {
	
	private final CruxGui cruxGui;
	private final JButton spectraButton = new JButton("Choose Spectra File");
	private final JTextField spectraFilename = new JTextField();
	
	private static Logger logger = 
		Logger.getLogger("edu.washington.gs.noble.crux.gui");
	
	class SpectraButtonListener implements ActionListener {
		public void actionPerformed(final ActionEvent event) {
			CruxAnalysisModel model = cruxGui.getAnalysisModel();
			final JFileChooser chooser = new JFileChooser();
			chooser.setDialogTitle("Choose spectra file");
			String fileName = model.getSpectraFilename();
			if (fileName != null && fileName.length() > 0) {
				chooser.setCurrentDirectory(new File(fileName));
			}
			else {
				chooser.setCurrentDirectory(new File("."));
			}
			final int result = chooser.showDialog(getParent(), "Choose spectra file");
			if (result == JFileChooser.APPROVE_OPTION){
				spectraFilename.setText(chooser.getSelectedFile().getPath());
				logger.info("User selected MS2 spectra file " + spectraFilename.getText() + ".");
			}
		}
	}
	
	public SpectraPanel(CruxGui cruxGui) {
		this.cruxGui = cruxGui;
		setLayout(new BoxLayout(this, BoxLayout.LINE_AXIS));
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
	
	public void loadDefaults() {
		CruxAnalysisModel model = cruxGui.getAnalysisModel();
		spectraFilename.setText(model.getDefaultSpectraFilename());
	}
	
	public void saveToModel() {
		CruxAnalysisModel model = cruxGui.getAnalysisModel();
		model.setSpectraFilename(spectraFilename.getText());
	}
	
	public void updateFromModel() {
		CruxAnalysisModel model = cruxGui.getAnalysisModel();
		spectraFilename.setText(model.getSpectraFilename());
	}
}