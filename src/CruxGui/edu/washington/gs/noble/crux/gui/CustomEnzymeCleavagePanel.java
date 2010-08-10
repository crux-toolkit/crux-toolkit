package edu.washington.gs.noble.crux.gui;

import java.util.HashMap;
import java.util.HashSet;
import java.awt.Color;
import java.awt.Component;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.Dimension;
import java.util.logging.Logger;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.JRadioButton;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JCheckBox;
import javax.swing.BorderFactory;

/**
 * The CustomEnzymeCleavagePanel class is the GUI elements of creating custom
 * enzyme cleavages. It allows users to specifiy customize which enzyme, for before and after,
 * and either to prevent or require.
 *
 * @author Michael Mathews
 *
 */


@SuppressWarnings("serial")
public class CustomEnzymeCleavagePanel extends CruxAdvancedParameterControl {
	private static Logger logger = 
		Logger.getLogger("edu.washington.gs.noble.crux.gui");

	private final CruxGui cruxGui;
	private final JLabel beforeClevageLabel = new JLabel("Before Cleavage:");
	private final JLabel afterClevageLabel = new JLabel("After Cleavage:");
	private HashMap<CruxAnalysisModel.AminoAcid, JCheckBox> beforeCheckBoxes = 
		new HashMap<CruxAnalysisModel.AminoAcid, JCheckBox>(CruxAnalysisModel.AminoAcid.getSize());
	private HashMap<CruxAnalysisModel.AminoAcid, JCheckBox> afterCheckBoxes =
		new HashMap<CruxAnalysisModel.AminoAcid, JCheckBox>(CruxAnalysisModel.AminoAcid.getSize());
	private final JCheckBox enableCustomEnzyme = new JCheckBox("Use Custom Enzyme");
	private final JLabel customEnzymeLabel = new JLabel("Custom Enzyme:");
	private final ButtonGroup beforeEnzymeTypeButtons = new ButtonGroup();
	private final JRadioButton beforeRequiredButton = new JRadioButton("Required");
	private final JRadioButton beforePreventButton = new JRadioButton("Prevent");
	private final ButtonGroup afterEnzymeTypeButtons = new ButtonGroup();
	private final JRadioButton afterRequiredButton = new JRadioButton("Required");
	private final JRadioButton afterPreventButton = new JRadioButton("Prevent");
	private final JPanel mainContainer = new JPanel();
	
	public CustomEnzymeCleavagePanel(CruxGui cruxGui) {
		super();
		this.cruxGui = cruxGui;
		this.setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
		this.setAlignmentX(Component.LEFT_ALIGNMENT);
		
		mainContainer.setBorder(BorderFactory.createLineBorder(Color.black));
		mainContainer.setBackground(Color.white);
		mainContainer.setLayout(new BoxLayout(mainContainer, BoxLayout.Y_AXIS));
		//mainContainer.setAlignmentX(Component.LEFT_ALIGNMENT);
		
		enableCustomEnzyme.addActionListener(new CustomEnzymeRadioListener());
		JPanel container = new JPanel();
		container.setLayout(new BoxLayout(container, BoxLayout.X_AXIS));
		container.setAlignmentY(Component.LEFT_ALIGNMENT);
		container.add(customEnzymeLabel);
		container.add(Box.createRigidArea(new Dimension(10,4)));
		container.add(enableCustomEnzyme);
		add(container);
		add(mainContainer);
		
		addComponents(beforeCheckBoxes, beforeClevageLabel, beforeRequiredButton, beforePreventButton, beforeEnzymeTypeButtons);
		addComponents(afterCheckBoxes, afterClevageLabel, afterRequiredButton, afterPreventButton, afterEnzymeTypeButtons);
	}
	
	private void addComponents(HashMap<CruxAnalysisModel.AminoAcid, JCheckBox> aminoAcidCheckBoxes, 
			JLabel enzymeLabel,
			JRadioButton requiredButton,
			JRadioButton preventButton,
			ButtonGroup enzymeTypeButtons){
		/* Label */
		JPanel container = getNewContainer();
		container.add(enzymeLabel);
		enzymeLabel.setAlignmentX(Component.LEFT_ALIGNMENT);
		enzymeLabel.setBackground(Color.white);
		mainContainer.add(container);
		/* Radio Buttons */
		container = getNewContainer();
		enzymeTypeButtons.add(requiredButton);
		enzymeTypeButtons.add(preventButton);
		requiredButton.setEnabled(false);
		preventButton.setEnabled(false);
		requiredButton.setSelected(true);
		container.add(requiredButton);
		container.add(preventButton);
		requiredButton.setBackground(Color.white);
		preventButton.setBackground(Color.white);
		requiredButton.setAlignmentX(Component.LEFT_ALIGNMENT);
		preventButton.setAlignmentX(Component.LEFT_ALIGNMENT);
		mainContainer.add(container);
		/* Check Boxes*/
		mainContainer.add(Box.createRigidArea(new Dimension(10,4)));
		container = getNewContainer();
		int i = 0;
		for (CruxAnalysisModel.AminoAcid aminoAcid : CruxAnalysisModel.AminoAcid.values()){	
			JCheckBox aminoAcidCheckBox = new JCheckBox(aminoAcid.toString());
			aminoAcidCheckBoxes.put(aminoAcid , aminoAcidCheckBox);
			container.add(aminoAcidCheckBox);
			container.add(Box.createRigidArea(new Dimension(10,4)));
			aminoAcidCheckBox.setBackground(Color.white);
			aminoAcidCheckBox.setEnabled(false);
			if (i % 5 == 4){
				mainContainer.add(container);
				container = getNewContainer();
			}
			i++;
		}
	}

	public void updateFromModel() {
		CruxAnalysisModel model = cruxGui.getAnalysisModel();
		beforePreventButton.setSelected(model.getPreventCustomEnzymeBefore());
		afterPreventButton.setSelected(model.getPreventCustomEnzymeAfter());
		enableCustomEnzyme.setSelected(model.getEnableCustomEnzyme());
		HashSet<CruxAnalysisModel.AminoAcid> selectedBefore = model.getSelectedCustomEnzymeBefore();
		HashSet<CruxAnalysisModel.AminoAcid> selectedAfter = model.getSelectedCustomEnzymeAfter();
		for (CruxAnalysisModel.AminoAcid aminoAcid : CruxAnalysisModel.AminoAcid.values()){
		if (selectedBefore.contains(aminoAcid)){
			beforeCheckBoxes.get(aminoAcid).setSelected(true);
			logger.info("Read the amino acid: "+aminoAcid.toString()+" was selected in the model for before cleavage");
		} else {
			beforeCheckBoxes.get(aminoAcid).setSelected(false);
		}
			if (selectedAfter.contains(aminoAcid)){
				afterCheckBoxes.get(aminoAcid).setSelected(true);
				logger.info("Read the amino acid: "+aminoAcid.toString()+" was selected in the model for after cleavage");
			} else {
				afterCheckBoxes.get(aminoAcid).setSelected(false);
			}
		}
		disableAminoSelection();	
		
	}
	
	public void loadDefaults() {
		CruxAnalysisModel model = cruxGui.getAnalysisModel();
		beforePreventButton.setSelected(model.getDefaultPreventCustomEnzymeBefore());
		afterPreventButton.setSelected(model.getDefaultPreventCustomEnzymeAfter());
		enableCustomEnzyme.setSelected(model.getDefaultEnableCustomEnzyme());
		HashSet<CruxAnalysisModel.AminoAcid> selectedBefore = model.getDefaultSelectedCustomEnzymeBefore();
		HashSet<CruxAnalysisModel.AminoAcid> selectedAfter = model.getDefaultSelectedCustomEnzymeAfter();
		for (CruxAnalysisModel.AminoAcid aminoAcid : CruxAnalysisModel.AminoAcid.values()){
			if (selectedBefore.contains(aminoAcid)){
				beforeCheckBoxes.get(aminoAcid).setSelected(true);
				logger.info("Read the amino acid: "+aminoAcid.toString()+" was selected in the model for before cleavage");
			} else {
				beforeCheckBoxes.get(aminoAcid).setSelected(false);
			}
			if (selectedAfter.contains(aminoAcid)){
				afterCheckBoxes.get(aminoAcid).setSelected(true);
				logger.info("Read the amino acid: "+aminoAcid.toString()+" was selected in the model for after cleavage");
			} else {
				afterCheckBoxes.get(aminoAcid).setSelected(false);
			}
		}
		disableAminoSelection();
	}
	
	public void saveToModel() {
		CruxAnalysisModel model = cruxGui.getAnalysisModel();
		HashSet<CruxAnalysisModel.AminoAcid> selectedBefore = 
				new HashSet<CruxAnalysisModel.AminoAcid>(CruxAnalysisModel.AminoAcid.getSize());
		HashSet<CruxAnalysisModel.AminoAcid> selectedAfter = 
				new HashSet<CruxAnalysisModel.AminoAcid>(CruxAnalysisModel.AminoAcid.getSize());
		for (CruxAnalysisModel.AminoAcid aminoAcid : CruxAnalysisModel.AminoAcid.values()){
			if (beforeCheckBoxes.get(aminoAcid).isSelected()){
				selectedBefore.add(aminoAcid);
				logger.info("Saving that the amino acid "+aminoAcid.toString()+" was selected for before cleavage to the model");
			} 
			if (afterCheckBoxes.get(aminoAcid).isSelected()){
				selectedAfter.add(aminoAcid);
				logger.info("Saving that the amino acid "+aminoAcid.toString()+" was selected for after cleavage to the model");
			} 
		}
		boolean preventBefore = beforePreventButton.isSelected();
		boolean preventAfter = afterPreventButton.isSelected();
		boolean enable = enableCustomEnzyme.isSelected();
		model.setCustomEnzymeCleavage(selectedBefore, selectedAfter, preventBefore, preventAfter, enable);
	}
	
	/* 
	 * Checks to see if user wants to use custom enzyme in gui and disables or enables 
	 * based on that value
	 * 
	 * */
	private void disableAminoSelection(){
		boolean enable = enableCustomEnzyme.isSelected();
		for (CruxAnalysisModel.AminoAcid aminoAcid : CruxAnalysisModel.AminoAcid.values()){
			beforeCheckBoxes.get(aminoAcid).setEnabled(enable);
			afterCheckBoxes.get(aminoAcid).setEnabled(enable);
			afterRequiredButton.setEnabled(enable);
			afterPreventButton.setEnabled(enable);
			beforeRequiredButton.setEnabled(enable);
			beforePreventButton.setEnabled(enable);
		}
	}
	
	private JPanel getNewContainer(){
		JPanel container = new JPanel();
		container.setLayout(new BoxLayout(container, BoxLayout.X_AXIS));
		container.setBackground(Color.white);
		container.setAlignmentY(Component.RIGHT_ALIGNMENT);
		return container;
	}
	
	class CustomEnzymeRadioListener implements ActionListener {
		public void actionPerformed(ActionEvent e) {
			disableAminoSelection();
		}
	}
}

