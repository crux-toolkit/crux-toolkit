package edu.washington.gs.noble.crux.gui;

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.util.logging.Logger;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JPanel;
import javax.swing.JCheckBox;


/**
 * The MissedCleavagesPanel is GUI element for user to 
 * choose whether they want to allow missed cleavages
 *
 * @author Michael Mathews
 *
 */
@SuppressWarnings("serial")
class MissedCleavagesPanel extends  CruxParameterControl {
	
	private static Logger logger = 
		Logger.getLogger("edu.washington.gs.noble.crux.gui");
	private final CruxGui cruxGui;
	private final JCheckBox allowMissedCleavages = new JCheckBox("Allow Missed Cleavages");

	
	public MissedCleavagesPanel(CruxGui cruxGui) {
		super();
		this.cruxGui = cruxGui;
		this.setAlignmentX(Component.LEFT_ALIGNMENT);
		this.setLayout(new BoxLayout(this, BoxLayout.X_AXIS));
		allowMissedCleavages.setSelected(false);
		add(allowMissedCleavages);
	}
	
	public void loadDefaults() {
		CruxAnalysisModel model = cruxGui.getAnalysisModel();
		allowMissedCleavages.setSelected(model.getDefaultAllowMissedCleavages());
		logger.info("Missed Cleavages read from model: " + model.getDefaultAllowMissedCleavages());
	}
	
	public void saveToModel() {
		CruxAnalysisModel model = cruxGui.getAnalysisModel();
		model.setAllowMissedCleavages(getAllowMissedCleavages());
	}
	
	public void updateFromModel() {
		CruxAnalysisModel model = cruxGui.getAnalysisModel();
		allowMissedCleavages.setSelected(model.getAllowMissedCleavages());
		logger.info("Missed Cleavages read from  model: " + model.getAllowMissedCleavages());
	}
	
	public boolean getAllowMissedCleavages() {
	       return allowMissedCleavages.isSelected();
	}
}