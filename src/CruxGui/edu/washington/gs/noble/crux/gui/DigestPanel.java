package edu.washington.gs.noble.crux.gui;

import java.awt.Color;
import java.awt.Component;

import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;

@SuppressWarnings("serial")
class DigestPanel extends JPanel {

	private final ButtonGroup digestButtons = new ButtonGroup();
	private final JLabel digestLabel = new JLabel("Digestion");
	private final JRadioButton fullDigest = new JRadioButton("Full Digest");
	private final JRadioButton partialDigest = new JRadioButton("Partial Digest");

	public DigestPanel(CruxModel model) {
		fullDigest.setSelected(true);
		partialDigest.setSelected(false);
		digestButtons.add(fullDigest);
		digestButtons.add(partialDigest);
		setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
		setBorder(BorderFactory.createLineBorder(Color.black));
		setBackground(Color.white);
		add(digestLabel);
		fullDigest.setBackground(Color.white);
		add(fullDigest);
		partialDigest.setBackground(Color.white);
		add(partialDigest);
		setAlignmentX(Component.LEFT_ALIGNMENT);
	}	
}