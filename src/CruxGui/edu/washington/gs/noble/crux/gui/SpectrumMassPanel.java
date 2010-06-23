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
 * The SpectrumMassPanel class is the GUI elements used to get and set the 
 * min and max masses for the spectra.
 *
 * @author Charles E. Grant
 *
 */
public class SpectrumMassPanel extends CruxAdvancedParameterControl {

	private static Logger logger = 
		Logger.getLogger("edu.washington.gs.noble.crux.gui");

	private final CruxGui cruxGui;
	private final JLabel maxLabel = new JLabel("Max. spectrum mass:");
	private final JPanel maxPanel = new JPanel();
	private final JPanel minPanel = new JPanel();
	private final JSpinner maxSpectralMass;
	private final SpinnerNumberModel maxModel;
	private final JLabel minLabel = new JLabel("Min. spectrum mass:");
	private final JSpinner minSpectralMass;
	private final SpinnerNumberModel minModel;

	public SpectrumMassPanel(CruxGui cruxGui) {
		this.cruxGui = cruxGui;
		CruxAnalysisModel model = cruxGui.getAnalysisModel();
		maxModel = new SpinnerNumberModel(model.getSpectrumMaxMassDefault(), 0.0, 1000000000.0, 0.1);
		maxSpectralMass = new JSpinner(maxModel);
		minModel = new SpinnerNumberModel(model.getSpectrumMinMassDefault(), 0.0, 1000000000.0, 0.1);
		minSpectralMass = new JSpinner(minModel);
		setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
		setAlignmentX(Component.LEFT_ALIGNMENT);
		initMaxPanel();
		initMinPanel();
		updateFromModel();
	}	
	
	private void initMaxPanel() {
		maxLabel.setMaximumSize(new Dimension(150, 30));
		maxSpectralMass.setMaximumSize(new Dimension(170, 30));
		maxPanel.setLayout(new BoxLayout(maxPanel, BoxLayout.X_AXIS));
		maxPanel.add(maxLabel);
		maxPanel.add(Box.createRigidArea(new Dimension(0,6)));
		maxPanel.add(maxSpectralMass);
		add(maxPanel);
	}
	
	private void initMinPanel() {
		minLabel.setMaximumSize(new Dimension(150, 30));
		minSpectralMass.setMaximumSize(new Dimension(170, 30));
		minPanel.setLayout(new BoxLayout(minPanel, BoxLayout.X_AXIS));
		minPanel.add(minLabel);
		minPanel.add(Box.createRigidArea(new Dimension(0,6)));
		minPanel.add(minSpectralMass);
		add(minPanel);
	}
	
	public void loadDefaults() {
		CruxAnalysisModel model = cruxGui.getAnalysisModel();
		minSpectralMass.setValue(model.getMinMassDefault());
		maxSpectralMass.setValue(model.getMaxMassDefault());
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
		maxSpectralMass.setValue(cruxGui.getAnalysisModel().getMaxMass());
		logger.info("Loaded parameter 'max-mass' from model: " + getMaxMass());
		minSpectralMass.setValue(cruxGui.getAnalysisModel().getMinMass());
		logger.info("Loaded parameter 'min-mass' from model: " + getMinMass());
	}
	
	public double getMaxMass() {
		return ((Double) maxSpectralMass.getValue()).doubleValue();
	}
	
	public double getMinMass() {
		return ((Double) minSpectralMass.getValue()).doubleValue();
	}

}
