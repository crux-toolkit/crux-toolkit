package edu.washington.gs.noble.crux.gui;

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.io.File;
import java.util.logging.Logger;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JPanel;
import javax.swing.JTextField;

@SuppressWarnings("serial")
class OutputDirPanel extends CruxParameterControl {
	
	private final CruxGui cruxGui;
	private final JButton outputDirButton = new JButton("Choose Output Directory");
	private final JTextField outputDirName = new JTextField();
	
	private static Logger logger = 
		Logger.getLogger("edu.washington.gs.noble.crux.gui");
	
	public OutputDirPanel(final CruxGui cruxGui) {
		this.cruxGui = cruxGui;
		setLayout(new BoxLayout(this, BoxLayout.LINE_AXIS));
		setAlignmentX(Component.LEFT_ALIGNMENT);
		outputDirButton.setAlignmentX(Component.LEFT_ALIGNMENT);
		outputDirButton.setAlignmentY(Component.CENTER_ALIGNMENT);
		outputDirButton.addActionListener(new OutputDirButtonListener());
		add(outputDirButton);
		add(Box.createRigidArea(new Dimension(12,0)));
		outputDirName.setAlignmentX(Component.LEFT_ALIGNMENT);
		outputDirName.setAlignmentY(Component.CENTER_ALIGNMENT);
		add(outputDirName);
		outputDirName.setMaximumSize(new Dimension(200,25));
	}
	
	public String getOutputDir() {
		return outputDirName.getText();
	}

	public void setOutputDirName(String dirName) {
		outputDirName.setText(dirName);
	}
	
	public void loadDefaults() {
		CruxAnalysisModel model = cruxGui.getAnalysisModel();
		outputDirName.setText(model.getOutputDirDefault());
	}

	public void saveToModel() {
		CruxAnalysisModel model = cruxGui.getAnalysisModel();
		model.setOutputDir(outputDirName.getText());
	}
	
	public void updateFromModel() {
		CruxAnalysisModel model = cruxGui.getAnalysisModel();
		outputDirName.setText(model.getOutputDir());
	}

	class OutputDirButtonListener implements ActionListener {
		public void actionPerformed(final ActionEvent event) {
			final JFileChooser chooser = new JFileChooser();
			chooser.setDialogTitle("Choose output directory");
			chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
			chooser.setAcceptAllFileFilterUsed(false);
			
			CruxAnalysisModel model = cruxGui.getAnalysisModel();
			String fileName = model.getOutputDir();
			if (fileName != null && fileName.length() > 0) {
				chooser.setCurrentDirectory(new File(fileName));
			}
			else {
				chooser.setCurrentDirectory(new File("."));
			}
			final int result = chooser.showDialog(getParent(), "Choose output directory");
			if (result == JFileChooser.APPROVE_OPTION){
				outputDirName.setText(chooser.getSelectedFile().getPath());
			}
		}
	}
}