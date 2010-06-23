package edu.washington.gs.noble.crux.gui;

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.util.logging.Logger;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;

/**
 * The MassPanel class is the GUI elements used to get and set the 
 * min and max masses for the peptides.
 *
 * @author Charles E. Grant
 *
 */
@SuppressWarnings("serial")
class MassPanel extends CruxAdvancedParameterControl {

	private static Logger logger = 
		Logger.getLogger("edu.washington.gs.noble.crux.gui");

	private final CruxGui cruxGui;
	private final JLabel maxLabel = new JLabel("Max. peptide mass:");
	private final JPanel maxPanel = new JPanel();
	private final JPanel minPanel = new JPanel();
	private final JSpinner maxMass;
	private final SpinnerNumberModel maxModel;
	private final JLabel minLabel = new JLabel("Min. peptide mass:");
	private final JSpinner minMass;
	private final SpinnerNumberModel minModel;

	public MassPanel(CruxGui cruxGui) {
		this.cruxGui = cruxGui;
		CruxAnalysisModel model = cruxGui.getAnalysisModel();
		maxModel = new SpinnerNumberModel(model.getMaxMassDefault(), 0.0, 10000.0, 0.1);
		maxMass = new JSpinner(maxModel);
		minModel = new SpinnerNumberModel(model.getMinMassDefault(), 0.0, 10000.0, 0.1);
		minMass = new JSpinner(minModel);
		setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
		setAlignmentX(Component.LEFT_ALIGNMENT);
		initMaxPanel();
		initMinPanel();
		updateFromModel();
	}	
	
	private void initMaxPanel() {
		maxLabel.setMaximumSize(new Dimension(130, 30));
		maxMass.setMaximumSize(new Dimension(100, 30));
		maxPanel.setLayout(new BoxLayout(maxPanel, BoxLayout.X_AXIS));
		maxPanel.add(maxLabel);
		maxPanel.add(Box.createRigidArea(new Dimension(0,6)));
		maxPanel.add(maxMass);
		add(maxPanel);
	}
	
	private void initMinPanel() {
		minLabel.setMaximumSize(new Dimension(130, 30));
		minMass.setMaximumSize(new Dimension(100, 30));
		minPanel.setLayout(new BoxLayout(minPanel, BoxLayout.X_AXIS));
		minPanel.add(minLabel);
		minPanel.add(Box.createRigidArea(new Dimension(0,6)));
		minPanel.add(minMass);
		add(minPanel);
	}
	
	public void loadDefaults() {
		CruxAnalysisModel model = cruxGui.getAnalysisModel();
		minMass.setValue(model.getMinMassDefault());
		maxMass.setValue(model.getMaxMassDefault());
	}
	
	public void saveToModel() {
		CruxAnalysisModel model = cruxGui.getAnalysisModel();
		model.setMaxMass(getMaxMass());
		logger.info("Saved parameter 'max-mass' to model: " + getMaxMass());
		model.setMinMass(getMinMass());
		logger.info("Saved parameter 'min-mass' to model: " + getMinMass());
	}
	
	public void updateFromModel() {
		CruxAnalysisModel model = cruxGui.getAnalysisModel();
		maxMass.setValue(cruxGui.getAnalysisModel().getMaxMass());
		logger.info("Loaded parameter 'max-mass' from model: " + getMaxMass());
		minMass.setValue(cruxGui.getAnalysisModel().getMinMass());
		logger.info("Loaded parameter 'min-mass' from model: " + getMinMass());
	}
	
	public double getMaxMass() {
		return ((Double) maxMass.getValue()).doubleValue();
	}
	
	public double getMinMass() {
		return ((Double) minMass.getValue()).doubleValue();
	}

}