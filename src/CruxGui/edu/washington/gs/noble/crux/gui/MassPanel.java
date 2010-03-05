package edu.washington.gs.noble.crux.gui;

import java.awt.Color;
import java.awt.Component;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.logging.Logger;

import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;

@SuppressWarnings("serial")
class MassPanel extends JPanel {

	private static Logger logger = 
		Logger.getLogger("edu.washington.gs.noble.crux.gui");

	private final ButtonGroup massButtons = new ButtonGroup();
	private final JLabel isotopicMassLabel = new JLabel("Isotopic Mass");
	private final JRadioButton averageMass = new JRadioButton("Average Mass");
	private final JRadioButton monoMass = new JRadioButton("Monoisotopic Mass");
	private CruxAnalysisModel model;

	public MassPanel(CruxAnalysisModel model) {
		this.model = model;
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
		averageMass.addItemListener(new IsotopicMassTypeChangeListener());
		monoMass.addItemListener(new IsotopicMassTypeChangeListener());
		updateFromModel();
	}	
	
	public void updateFromModel() {
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
	
	class IsotopicMassTypeChangeListener implements ItemListener {
		public void itemStateChanged(final ItemEvent event) {
			if (averageMass.isSelected()) {
				model.setMassType(CruxAnalysisModel.IsotopicMassType.AVERAGE);
			}
			else if (monoMass.isSelected()) {
				model.setMassType(CruxAnalysisModel.IsotopicMassType.MONOISOTOPIC);		
			}
		}
	}
}