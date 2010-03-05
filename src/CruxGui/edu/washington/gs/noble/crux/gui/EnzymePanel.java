package edu.washington.gs.noble.crux.gui;

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.logging.Logger;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;

@SuppressWarnings("serial")
class EnzymePanel extends JPanel {
	
	private static Logger logger = 
		Logger.getLogger("edu.washington.gs.noble.crux.gui");

	CruxAnalysisModel model = null;
	final JLabel enzymeLabel = new JLabel("Enzyme:");
	final JComboBox enzymeCombo = new JComboBox(CruxAnalysisModel.Enzyme.values());
	
	public EnzymePanel(CruxAnalysisModel model) {
		super();
		this.model = model;
		this.setBackground(Color.white);
		this.setAlignmentX(Component.LEFT_ALIGNMENT);
		this.setLayout(new BoxLayout(this, BoxLayout.X_AXIS));
		enzymeLabel.setBackground(Color.white);
		add(enzymeLabel);
		add(Box.createRigidArea(new Dimension(6,0)));
		enzymeCombo.setSelectedIndex(0);
		enzymeCombo.setBackground(Color.white);
		enzymeCombo.setMaximumSize(new Dimension(236, 24));
		add(enzymeCombo);
		enzymeCombo.addItemListener(new EnzymeChangeListener());
	}
	
	void updateFromModel() {
		enzymeCombo.setSelectedIndex(model.getEnzyme().ordinal());
		logger.info("Enzyme read from model: " + model.getEnzyme().toString());
	}
	
	class EnzymeChangeListener implements ItemListener {
		public void itemStateChanged(final ItemEvent event) {
			CruxAnalysisModel.Enzyme e = (CruxAnalysisModel.Enzyme) enzymeCombo.getSelectedItem();
			model.setEnzyme(e);
		}
	}	
}