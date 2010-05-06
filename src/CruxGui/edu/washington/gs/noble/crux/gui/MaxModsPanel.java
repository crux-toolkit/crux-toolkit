package edu.washington.gs.noble.crux.gui;

import java.awt.Component;
import java.awt.Dimension;
import java.util.logging.Logger;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;

/**
 * The MaxModsPanel class is the GUI elements used to get and set the 
 * max number of modfications for the peptides.
 *
 * @author Charles E. Grant
 *
 */
@SuppressWarnings("serial")
public class MaxModsPanel extends CruxParameterControl {
		private static Logger logger = 
			Logger.getLogger("edu.washington.gs.noble.crux.gui");

		private final CruxGui cruxGui;
		private final JLabel maxLabel = new JLabel("Max. peptide mods:");
		private final JPanel maxPanel = new JPanel();
		private final JSpinner maxMods;
		private final SpinnerNumberModel maxModel;

		public MaxModsPanel(CruxGui cruxGui) {
			this.cruxGui = cruxGui;
			CruxAnalysisModel model = cruxGui.getAnalysisModel();
			maxModel = new SpinnerNumberModel(model.getMaxModsDefault(), 0, 255, 1);
			maxMods = new JSpinner(maxModel);
			setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
			setAlignmentX(Component.LEFT_ALIGNMENT);
			initMaxPanel();
			updateFromModel();
		}	
		
		private void initMaxPanel() {
			maxLabel.setMaximumSize(new Dimension(130, 30));
			maxMods.setMaximumSize(new Dimension(130, 30));
			maxPanel.setLayout(new BoxLayout(maxPanel, BoxLayout.X_AXIS));
			maxPanel.add(maxLabel);
			maxPanel.add(Box.createRigidArea(new Dimension(0,6)));
			maxPanel.add(maxMods);
			add(maxPanel);
		}
		
		public void loadDefaults() {
			CruxAnalysisModel model = cruxGui.getAnalysisModel();
			maxMods.setValue(model.getMaxModsDefault());
		}
		
		public void saveToModel() {
			CruxAnalysisModel model = cruxGui.getAnalysisModel();
			model.setMaxMods(getMaxMods());
			logger.info("Saved parameter 'max-mods' to model: " + getMaxMods());
		}
		
		private int getMaxMods() {
			return ((Integer) maxMods.getValue()).intValue();
		}

		public void updateFromModel() {
			CruxAnalysisModel model = cruxGui.getAnalysisModel();
			maxMods.setValue(cruxGui.getAnalysisModel().getMaxMods());
			logger.info("Loaded parameter 'max-mods' from model: " + getMaxMods());
		}
		

}
