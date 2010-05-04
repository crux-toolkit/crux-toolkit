package edu.washington.gs.noble.crux.gui;

import java.awt.Color;
import java.awt.Component;
import java.awt.event.ItemListener;
import java.util.logging.Logger;

import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;

/**
 * The MassTypePanel class is the GUI elements used to get and set the 
 * mass type for the peptides (average or monoisotopic).
 *
 * @author Charles E. Grant
 *
 */
@SuppressWarnings("serial")
class MassTypePanel extends CruxParameterControl {

	private static Logger logger = 
		Logger.getLogger("edu.washington.gs.noble.crux.gui");

	private final CruxGui cruxGui;
	private final ButtonGroup massButtons = new ButtonGroup();
	private final JLabel isotopicMassLabel = new JLabel("Isotopic Mass");
	private final JRadioButton averageMass = new JRadioButton("Average Mass");
	private final JRadioButton monoMass = new JRadioButton("Monoisotopic Mass");

	public MassTypePanel(CruxGui cruxGui) {
		this.cruxGui = cruxGui;
		massButtons.add(averageMass);
		massButtons.add(monoMass);
		setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
		setBorder(BorderFactory.createLineBorder(Color.black));
		setBackground(Color.white);
		add(isotopicMassLabel);
		averageMass.setBackground(Color.white);
		add(averageMass);
		monoMass.setBackground(Color.white);
		add(monoMass);
		setAlignmentX(Component.LEFT_ALIGNMENT);
		updateFromModel();
	}	
	
	public void loadDefaults() {
		CruxAnalysisModel model = cruxGui.getAnalysisModel();
		if (model.getDefaultMassType() == CruxAnalysisModel.IsotopicMassType.AVERAGE) {
			averageMass.setSelected(true);
			monoMass.setSelected(false);
			logger.info("Read parameter 'IsotopicMassType' from model: AVERAGE");
		}
		else if (model.getDefaultMassType() == CruxAnalysisModel.IsotopicMassType.MONOISOTOPIC) {
			averageMass.setSelected(false);
			monoMass.setSelected(true);
			logger.info("Read default parameter 'IsotopicMassType' from model: MONOISOTOPIC");
		}
		else {
			logger.info("Read unrecognized value of parameter 'IsotopicMassType' from model: " + model.getDefaultMassType().toString() + ".");
		}
	}
	
	public void saveToModel() {
		CruxAnalysisModel model = cruxGui.getAnalysisModel();
		model.setMassType(getMassType());
	}
	
	public void updateFromModel() {
		CruxAnalysisModel model = cruxGui.getAnalysisModel();
		if (model.getMassType() == CruxAnalysisModel.IsotopicMassType.AVERAGE) {
			averageMass.setSelected(true);
			monoMass.setSelected(false);
			logger.info("Read parameter 'IsotopicMassType' from model: AVERAGE");
		}
		else if (model.getMassType() == CruxAnalysisModel.IsotopicMassType.MONOISOTOPIC) {
			averageMass.setSelected(false);
			monoMass.setSelected(true);
			logger.info("Read parameter 'IsotopicMassType' from model: MONOISOTOPIC");
		}
		else {
			logger.info("Read unrecognized value of parameter 'IsotopicMassType' from model: " + model.getMassType().toString() + ".");
		}
	}
	
	public CruxAnalysisModel.IsotopicMassType getMassType() {
		if (averageMass.isSelected()) {
			return CruxAnalysisModel.IsotopicMassType.AVERAGE;
		}
		else if (monoMass.isSelected()) {
			return CruxAnalysisModel.IsotopicMassType.MONOISOTOPIC;
		}
		else {
			return null;
		}
	}

}