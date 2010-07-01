package edu.washington.gs.noble.crux.gui;

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.util.logging.Logger;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;



/**
 * The EnzymePanel class is the GUI elements used to get and set the 
 * enzyme used to perform the digestion of the peptides.
 *
 * @author Charles E. Grant
 *
 */
@SuppressWarnings("serial")
class EnzymePanel extends  CruxAdvancedParameterControl {
	
	private static Logger logger = 
		Logger.getLogger("edu.washington.gs.noble.crux.gui");

	private final CruxGui cruxGui;
	private final JLabel enzymeLabel = new JLabel("Enzyme:");
	private final JComboBox enzymeCombo = new JComboBox(CruxAnalysisModel.Enzyme.values());
	
	public EnzymePanel(CruxGui cruxGui) {
		super();
		this.cruxGui = cruxGui;
		this.setAlignmentX(Component.LEFT_ALIGNMENT);
		this.setLayout(new BoxLayout(this, BoxLayout.X_AXIS));
		add(enzymeLabel);
		add(Box.createRigidArea(new Dimension(6,0)));
		enzymeCombo.setSelectedIndex(0);
		enzymeCombo.setBackground(Color.white);
		enzymeCombo.setMaximumSize(new Dimension(236, 24));
		add(enzymeCombo);
	}
	
	public void loadDefaults() {
		CruxAnalysisModel model = cruxGui.getAnalysisModel();
		enzymeCombo.setSelectedIndex(model.getDefaultEnzyme().ordinal());
		logger.info("Enzyme read from model: " + model.getDefaultEnzyme().toString());
	}
	
	public void saveToModel() {
		CruxAnalysisModel model = cruxGui.getAnalysisModel();
		model.setEnzyme(getEnzyme());
	}
	
	public void updateFromModel() {
		CruxAnalysisModel model = cruxGui.getAnalysisModel();
		enzymeCombo.setSelectedIndex(model.getEnzyme().ordinal());
		logger.info("Enzyme read from model: " + model.getEnzyme().toString());
	}
	
	public CruxAnalysisModel.Enzyme getEnzyme() {
		return (CruxAnalysisModel.Enzyme) enzymeCombo.getSelectedItem();
	}
}