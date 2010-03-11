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
class ProteinDBPanel extends JPanel {
	
	CruxAnalysisModel model;
	JButton proteinDBButton = new JButton("Choose Protein Database");
	JTextField proteinFileName = new JTextField();
	
	private static Logger logger = 
		Logger.getLogger("edu.washington.gs.noble.crux.gui");
	
	public ProteinDBPanel(CruxAnalysisModel model) {
		this.model = model;
		setLayout(new BoxLayout(this, BoxLayout.LINE_AXIS));
		setBackground(Color.white);
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
	
	public void updateFromModel() {
		proteinFileName.setText(model.getProteinDatabase());
	}

	class ProteinDBButtonListener implements ActionListener {
		public void actionPerformed(final ActionEvent event) {
			final JFileChooser chooser = new JFileChooser();
			String fileName = model.getProteinDatabase();
			if (fileName != null && fileName.length() > 0) {
				chooser.setCurrentDirectory(new File(fileName));
			}
			else {
				chooser.setCurrentDirectory(new File("."));
			}
			final int result = chooser.showOpenDialog(getParent());
			if (result == JFileChooser.APPROVE_OPTION){
				proteinFileName.setText(chooser.getSelectedFile().getPath());
			}
		}
	}
}