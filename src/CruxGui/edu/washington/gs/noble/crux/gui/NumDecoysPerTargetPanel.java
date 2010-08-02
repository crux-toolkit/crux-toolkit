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
 * The DecoysPerTarget class is the GUI elements used to get and set the 
 * number of decoys peptides per target peptide for searching for matches.
 *
 * @author Michael Mathews
 *
 */
@SuppressWarnings("serial")
class NumDecoysPerTargetPanel extends CruxParameterControl {

	private static Logger logger = 
		Logger.getLogger("edu.washington.gs.noble.crux.gui");

	private final CruxGui cruxGui;
	private final JLabel numDecoysPerTargetLabel = new JLabel("Number of decoys per target:");
	private final JPanel numDecoysPerTargetPanel = new JPanel();
	private final JSpinner numDecoysPerTarget;
	private final SpinnerNumberModel numDecoysPerTargetModel;

	public NumDecoysPerTargetPanel(CruxGui cruxGui) {
		this.cruxGui = cruxGui;
		CruxAnalysisModel model = cruxGui.getAnalysisModel();
		numDecoysPerTargetModel = new SpinnerNumberModel(model.getNumDecoysPerTargetDefault(), 1, 10000, 1);
		numDecoysPerTarget = new JSpinner(numDecoysPerTargetModel);
		setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
		setAlignmentX(Component.LEFT_ALIGNMENT);
		initDecoysPerTargetPanel();
		updateFromModel();
	}	
	
	private void initDecoysPerTargetPanel() {
		numDecoysPerTargetLabel.setMaximumSize(new Dimension(200, 30));
		numDecoysPerTarget.setMaximumSize(new Dimension(120, 30));
		numDecoysPerTargetPanel.setLayout(new BoxLayout(numDecoysPerTargetPanel, BoxLayout.X_AXIS));
		numDecoysPerTargetPanel.add(numDecoysPerTargetLabel);
		numDecoysPerTargetPanel.add(Box.createRigidArea(new Dimension(0,6)));
		numDecoysPerTargetPanel.add(numDecoysPerTarget);
		add(numDecoysPerTargetPanel);
	}
	
	
	public void loadDefaults() {
		CruxAnalysisModel model = cruxGui.getAnalysisModel();
		numDecoysPerTarget.setValue(model.getNumDecoysPerTargetDefault());
	}
	
	public void saveToModel() {
		CruxAnalysisModel model = cruxGui.getAnalysisModel();
		model.setNumDecoysPerTarget(getNumDecoysPerTarget());
		logger.info("Saved parameter 'numDecoysPerTarget' to model: " + model.getNumDecoysPerTarget());
	}
	
	public void updateFromModel() {
		CruxAnalysisModel model = cruxGui.getAnalysisModel();
		numDecoysPerTarget.setValue(model.getNumDecoysPerTarget());
		logger.info("Loaded parameter 'numDecoysPerTarget' from model: " + model.getNumDecoysPerTarget());
	}
	
	public int getNumDecoysPerTarget() {
	    return ((Integer) numDecoysPerTarget.getValue()).intValue();
	}
	

}