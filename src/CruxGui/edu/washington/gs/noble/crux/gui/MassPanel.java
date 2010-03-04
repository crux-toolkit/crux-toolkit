package edu.washington.gs.noble.crux.gui;

import java.awt.Color;
import java.awt.Component;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.logging.Logger;

import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;

@SuppressWarnings("serial")
class MassPanel extends JPanel implements ActionListener {

	private static Logger logger = 
		Logger.getLogger("edu.washington.gs.noble.crux.gui");

	private final ButtonGroup massButtons = new ButtonGroup();
	private final JLabel isotopicMassLabel = new JLabel("Isotopic Mass");
	private final JRadioButton averageMass = new JRadioButton("Average Mass");
	private final JRadioButton monoMass = new JRadioButton("Monoisotopic Mass");
	private CruxModel model;

	public MassPanel(CruxModel model) {
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
		updateFromModel();
		averageMass.addActionListener(this);
		monoMass.addActionListener(this);
	}	
	
	public void actionPerformed(ActionEvent event) {
		if (averageMass.isSelected()) {
			model.setMassType(CruxModel.IsotopicMassType.AVERAGE);
			logger.info("Model IsotopicMassType set to AVERAGE");
		}
		else if (monoMass.isSelected()) {
			model.setMassType(CruxModel.IsotopicMassType.MONOISOTOPIC);		
			logger.info("Model IsotopicMassType set to MONOISOTOPIC");
		}
		
	}
	
	public void updateFromModel() {
		if (model.getMassType() == CruxModel.IsotopicMassType.AVERAGE) {
			averageMass.setSelected(true);
			monoMass.setSelected(false);
			logger.info("Mass type set to AVERAGE");
		}
		else if (model.getMassType() == CruxModel.IsotopicMassType.MONOISOTOPIC) {
			averageMass.setSelected(false);
			monoMass.setSelected(true);
			logger.info("Mass type set to MONOISOTOPIC");
		}
		else {
			logger.info("Unknown mass type " + model.getMassType().toString());
		}
	}
}