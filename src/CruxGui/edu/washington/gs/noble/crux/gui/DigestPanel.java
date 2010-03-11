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
class DigestPanel extends JPanel {

	private static Logger logger = 
		Logger.getLogger("edu.washington.gs.noble.crux.gui");

	private final ButtonGroup digestButtons = new ButtonGroup();
	private final JLabel digestLabel = new JLabel("Digestion");
	private final JRadioButton fullDigest = new JRadioButton("Full Digest");
	private final JRadioButton partialDigest = new JRadioButton("Partial Digest");
	private CruxAnalysisModel model;

	public DigestPanel(CruxAnalysisModel model) {
		this.model = model;
		fullDigest.setSelected(true);
		partialDigest.setSelected(false);
		digestButtons.add(fullDigest);
		digestButtons.add(partialDigest);
		setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
		setBorder(BorderFactory.createLineBorder(Color.black));
		setBackground(Color.white);
		setAlignmentX(Component.LEFT_ALIGNMENT);
		add(digestLabel);
		fullDigest.setBackground(Color.white);
		add(fullDigest);
		partialDigest.setBackground(Color.white);
		add(partialDigest);
		updateFromModel();
	}	
	
	void updateFromModel() {
		if (model.getDigestType() == CruxAnalysisModel.DigestType.FULL) {
			fullDigest.setSelected(true);
			partialDigest.setSelected(false);
			logger.info("Read parameter 'DigestType' from model: " + CruxAnalysisModel.DigestType.FULL.toString());
		}
		else if (model.getDigestType() == CruxAnalysisModel.DigestType.PARTIAL) {
			fullDigest.setSelected(false);
			partialDigest.setSelected(true);
			logger.info("Read parameter 'DigestType' from model: " + CruxAnalysisModel.DigestType.PARTIAL.toString());
		}
		else {
			logger.info("Read unrecognized value of parameter 'DigestType' from model: " + model.getDigestType().toString());
		}

	}
	
	public CruxAnalysisModel.DigestType getDigestType() {
		if (fullDigest.isSelected()) {
			 return CruxAnalysisModel.DigestType.FULL;
		}
		else if (partialDigest.isSelected()) {
			return CruxAnalysisModel.DigestType.PARTIAL;		
		}
		else {
			return null;
		}
	}
}