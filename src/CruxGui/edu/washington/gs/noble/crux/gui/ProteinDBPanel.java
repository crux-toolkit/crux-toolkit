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
class ProteinDBPanel extends CruxParameterControl {
	
	private final CruxGui cruxGui;
	private final JButton proteinDBButton = new JButton("Choose Protein Database");
	private final JTextField proteinFileName = new JTextField();
	
	private static Logger logger = 
		Logger.getLogger("edu.washington.gs.noble.crux.gui");
	
	public ProteinDBPanel(final CruxGui cruxGui) {
		this.cruxGui = cruxGui;
		setLayout(new BoxLayout(this, BoxLayout.LINE_AXIS));
		setAlignmentX(Component.LEFT_ALIGNMENT);
		proteinDBButton.setAlignmentX(Component.LEFT_ALIGNMENT);
		proteinDBButton.setAlignmentY(Component.CENTER_ALIGNMENT);
		proteinDBButton.addActionListener(new ProteinDBButtonListener());
		add(proteinDBButton);
		add(Box.createRigidArea(new Dimension(12,0)));
		proteinFileName.setAlignmentX(Component.LEFT_ALIGNMENT);
		proteinFileName.setAlignmentY(Component.CENTER_ALIGNMENT);
		add(proteinFileName);
		proteinFileName.setMaximumSize(new Dimension(200,25));
	}
	
	public String getProteinDatabase() {
		return proteinFileName.getText();
	}

	public void setProteinDatabase(String fileName) {
		proteinFileName.setText(fileName);
	}
	
	public void loadDefaults() {
		CruxAnalysisModel model = cruxGui.getAnalysisModel();
		proteinFileName.setText(model.getDefaultProteinDatabase());
	}

	public void saveToModel() {
		CruxAnalysisModel model = cruxGui.getAnalysisModel();
		model.setProteinDatabase(proteinFileName.getText());
	}
	
	public void updateFromModel() {
		CruxAnalysisModel model = cruxGui.getAnalysisModel();
		proteinFileName.setText(model.getProteinDatabase());
	}

	class ProteinDBButtonListener implements ActionListener {
		public void actionPerformed(final ActionEvent event) {
			final JFileChooser chooser = new JFileChooser();
			chooser.setDialogTitle("Choose protein database");
			CruxAnalysisModel model = cruxGui.getAnalysisModel();
			String fileName = model.getProteinDatabase();
			if (fileName != null && fileName.length() > 0) {
				chooser.setCurrentDirectory(new File(fileName));
			}
			else {
				chooser.setCurrentDirectory(new File("."));
			}
			final int result = chooser.showDialog(getParent(), "Choose protein database");
			if (result == JFileChooser.APPROVE_OPTION){
				proteinFileName.setText(chooser.getSelectedFile().getPath());
			}
		}
	}
}