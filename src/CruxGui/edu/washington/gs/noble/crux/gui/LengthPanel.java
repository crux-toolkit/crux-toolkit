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
 * The LengthPanel class is the GUI elements used to get and set the 
 * min and max lengths for the peptides.
 *
 * @author Charles E. Grant
 *
 */
@SuppressWarnings("serial")
class LengthPanel extends CruxAdvancedParameterControl {

	private static Logger logger = 
		Logger.getLogger("edu.washington.gs.noble.crux.gui");

	private final CruxGui cruxGui;
	private final JLabel maxLabel = new JLabel("Max. peptide length:");
	private final JPanel maxPanel = new JPanel();
	private final JPanel minPanel = new JPanel();
	private final JSpinner maxLength;
	private final SpinnerNumberModel maxModel;
	private final JLabel minLabel = new JLabel("Min. peptide length:");
	private final JSpinner minLength;
	private final SpinnerNumberModel minModel;

	public LengthPanel(CruxGui cruxGui) {
		this.cruxGui = cruxGui;
		CruxAnalysisModel model = cruxGui.getAnalysisModel();
		maxModel = new SpinnerNumberModel();
		maxLength = new JSpinner(maxModel);
		minModel = new SpinnerNumberModel();
		minLength = new JSpinner(minModel);
		setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
		setAlignmentX(Component.LEFT_ALIGNMENT);
		initMaxPanel();
		initMinPanel();
		updateFromModel();
	}	
	
	private void initMaxPanel() {
		maxLabel.setMaximumSize(new Dimension(130, 30));
		maxLength.setMaximumSize(new Dimension(75, 30));
		maxPanel.setLayout(new BoxLayout(maxPanel, BoxLayout.X_AXIS));
		maxPanel.add(maxLabel);
		maxPanel.add(Box.createRigidArea(new Dimension(0,20)));
		maxPanel.add(maxLength);
		add(maxPanel);
		add(Box.createRigidArea(new Dimension(10,4)));
	}
	
	private void initMinPanel() {
		minLabel.setMaximumSize(new Dimension(130, 30));
		minLength.setMaximumSize(new Dimension(75, 30));
		minPanel.setLayout(new BoxLayout(minPanel, BoxLayout.X_AXIS));
		minPanel.add(minLabel);
		minPanel.add(Box.createRigidArea(new Dimension(0,20)));
		minPanel.add(minLength);
		add(minPanel);
	}
	
	public void loadDefaults() {
		CruxAnalysisModel model = cruxGui.getAnalysisModel();
		minLength.setValue(model.getMinLengthDefault());
		maxLength.setValue(model.getMaxLengthDefault());
	}
	
	public void saveToModel() {
		CruxAnalysisModel model = cruxGui.getAnalysisModel();
		model.setMaxLength(getMaxLength());
		logger.info("Saved parameter 'max-length' to model: " + getMaxLength());
		model.setMinLength(getMinLength());
		logger.info("Saved parameter 'min-length' to model: " + getMinLength());
	}
	
	public void updateFromModel() {
		CruxAnalysisModel model = cruxGui.getAnalysisModel();
		maxLength.setValue(cruxGui.getAnalysisModel().getMaxLength());
		logger.info("Loaded parameter 'max-length' from model: " + getMaxLength());
		minLength.setValue(cruxGui.getAnalysisModel().getMinLength());
		logger.info("Loaded parameter 'min-length' from model: " + getMinLength());
	}
	
	public int getMaxLength() {
		return ((Integer) maxLength.getValue()).intValue();
	}
	
	public int getMinLength() {
		return ((Integer) minLength.getValue()).intValue();
	}

}