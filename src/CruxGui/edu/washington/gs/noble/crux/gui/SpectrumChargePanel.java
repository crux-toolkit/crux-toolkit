package edu.washington.gs.noble.crux.gui;

import java.awt.Color;
import java.awt.Component;
import java.util.logging.Logger;

import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.JLabel;
import javax.swing.JRadioButton;

/**
 * The SpectrumChargePanel class is the GUI elements used to get and set the 
 * spectrum charge: 1, 2, 3, or all.
 *
 * @author Charles E. Grant
 *
 */
public class SpectrumChargePanel extends CruxParameterControl {
	
	private static Logger logger = 
		Logger.getLogger("edu.washington.gs.noble.crux.gui");

	private final CruxGui cruxGui;
	private final ButtonGroup chargeButtons = new ButtonGroup();
	private final JLabel chargeLabel = new JLabel("Spectrum charge");
	private final JRadioButton charge1 = new JRadioButton("1");
	private final JRadioButton charge2 = new JRadioButton("2");
	private final JRadioButton charge3 = new JRadioButton("3");
	private final JRadioButton chargeAll = new JRadioButton("all");

	public SpectrumChargePanel(CruxGui cruxGui) {
		this.cruxGui = cruxGui;
		charge1.setSelected(false);
		charge1.setBackground(Color.white);
		chargeButtons.add(charge1);
		charge2.setSelected(false);
		charge2.setBackground(Color.white);
		chargeButtons.add(charge2);
		charge3.setSelected(false);
		charge3.setBackground(Color.white);
		chargeButtons.add(charge3);
		chargeAll.setSelected(true);
		chargeAll.setBackground(Color.white);
		chargeButtons.add(chargeAll);
		setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
		setBorder(BorderFactory.createLineBorder(Color.black));
		setBackground(Color.white);
		setAlignmentX(Component.LEFT_ALIGNMENT);
		add(chargeLabel);
		add(charge1);
		add(charge2);
		add(charge3);
		add(chargeAll);
		updateFromModel();
	}	
	
	
	public void loadDefaults() {
		CruxAnalysisModel model = cruxGui.getAnalysisModel();
		switch(model.getSpectrumChargeDefault()) {
			case ONE:
				charge1.setSelected(true);
				charge2.setSelected(false);
				charge3.setSelected(false);
				chargeAll.setSelected(false);
				break;
			case TWO:
				charge1.setSelected(false);
				charge2.setSelected(true);
				charge3.setSelected(false);
				chargeAll.setSelected(false);
				break;
			case THREE:
				charge1.setSelected(false);
				charge2.setSelected(false);
				charge3.setSelected(true);
				chargeAll.setSelected(false);
				break;
			case ALL:
				charge1.setSelected(false);
				charge2.setSelected(false);
				charge3.setSelected(false);
				chargeAll.setSelected(true);
				break;
			default:
				break;
				
		}
		logger.info("Read parameter 'Spectrum charge' from model: " + model.getSpectrumChargeDefault().toString());
	}
	
	public void saveToModel() {
		CruxAnalysisModel model = cruxGui.getAnalysisModel();
		model.setSpectrumCharge(getSpectrumCharge());
	}
	
	public void updateFromModel() {
		CruxAnalysisModel model = cruxGui.getAnalysisModel();
		switch(model.getSpectrumCharge()) {
			case ONE:
				charge1.setSelected(true);
				charge2.setSelected(false);
				charge3.setSelected(false);
				chargeAll.setSelected(false);
				break;
			case TWO:
				charge1.setSelected(false);
				charge2.setSelected(true);
				charge3.setSelected(false);
				chargeAll.setSelected(false);
				break;
			case THREE:
				charge1.setSelected(false);
				charge2.setSelected(false);
				charge3.setSelected(true);
				chargeAll.setSelected(false);
				break;
			case ALL:
				charge1.setSelected(false);
				charge2.setSelected(false);
				charge3.setSelected(false);
				chargeAll.setSelected(true);
				break;
			default:
				break;
		}
		logger.info("Read parameter 'Spectrum charge' from model: " + model.getSpectrumCharge().toString());
	}
	
	public CruxAnalysisModel.AllowedSpectrumCharge getSpectrumCharge() {
		if (charge1.isSelected()) {
			 return CruxAnalysisModel.AllowedSpectrumCharge.ONE;
		}
		else if (charge2.isSelected()) {
			 return CruxAnalysisModel.AllowedSpectrumCharge.TWO;
		}
		else if (charge3.isSelected()) {
			 return CruxAnalysisModel.AllowedSpectrumCharge.THREE;
		}
		else if (chargeAll.isSelected()) {
			 return CruxAnalysisModel.AllowedSpectrumCharge.ALL;
		}
		else {
			return null;
		}
	}
}
