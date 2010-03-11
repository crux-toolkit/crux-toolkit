package edu.washington.gs.noble.crux.gui;

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.logging.Logger;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JToggleButton;
import javax.swing.JCheckBox;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.border.EtchedBorder;

@SuppressWarnings("serial")
class IndexPanel extends JPanel implements ItemListener {
	
	private static Logger logger = 
		Logger.getLogger("edu.washington.gs.noble.crux.gui");
	
	private final JCheckBox runToolCheckBox = new JCheckBox("Run this tool");
	private final JCheckBox showAdvancedParameters = new JCheckBox("Show advanced parameters");
	private ProteinDBPanel proteinDBPanel = null;
 	private EnzymePanel enzymePanel = null;
	private MassPanel massPanel = null;
	private DigestPanel digestPanel = null;
	private final JCheckBox missedCleavages = new JCheckBox("Allow missed-cleavages");
	private final JButton saveButton = new JButton("Save");
	private final JButton cancelButton = new JButton("Cancel");
	private final JButton loadDefaults = new JButton("Load Defaults");
	private CruxAnalysisModel model = null;
	private CruxComponentButton button;
	private JToggleButton dummyButton;
	
	public IndexPanel(CruxAnalysisModel model, final CruxComponentButton button, final JToggleButton dummy) {
		super();
		this.model = model;
		this.button = button;
		this.dummyButton = dummy;
		this.setVisible(false);
		setBorder(
			BorderFactory.createTitledBorder(
				BorderFactory.createEtchedBorder(EtchedBorder.RAISED), 
				"create-index parameters"
			)
		);
		setBackground(Color.white);
		setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
		proteinDBPanel = new ProteinDBPanel(model);
		enzymePanel = new EnzymePanel(model);
		massPanel = new MassPanel(model);
		digestPanel = new DigestPanel(model);
		runToolCheckBox.setBackground(Color.white);
		missedCleavages.setBackground(Color.white);
		showAdvancedParameters.setBackground(Color.white);
		showAdvancedParameters.setEnabled(false);
		JPanel buttonPanel = new JPanel();
		buttonPanel.setLayout(new BoxLayout(buttonPanel, BoxLayout.LINE_AXIS));
		buttonPanel.setBackground(Color.white);
		buttonPanel.setAlignmentX(Component.LEFT_ALIGNMENT);
		saveButton.addActionListener(new SaveButtonListener());
		// saveButton.setEnabled(false);
		buttonPanel.add(Box.createRigidArea(new Dimension(36,0)));
		buttonPanel.add(saveButton);
		cancelButton.addActionListener(new CancelButtonListener());
		buttonPanel.add(Box.createRigidArea(new Dimension(12,0)));
		buttonPanel.add(cancelButton);
		buttonPanel.add(Box.createRigidArea(new Dimension(12,0)));
		loadDefaults.setEnabled(false);
		buttonPanel.add(loadDefaults);
		add(Box.createRigidArea(new Dimension(0,12)));
		add(runToolCheckBox);
		add(Box.createRigidArea(new Dimension(0,12)));
		add(showAdvancedParameters);
		add(Box.createRigidArea(new Dimension(0,12)));
	    add(proteinDBPanel);
		add(Box.createRigidArea(new Dimension(0,12)));
		add(massPanel);
		add(Box.createRigidArea(new Dimension(0,12)));
		add(digestPanel);
		add(Box.createRigidArea(new Dimension(0,12)));
		add(enzymePanel);
		add(Box.createRigidArea(new Dimension(0,12)));
		add(missedCleavages);
		add(Box.createRigidArea(new Dimension(0,12)));
		add(buttonPanel);
		runToolCheckBox.addItemListener(new RunComponentChangeListener());
		updateFromModel();
		setVisible(false);
	}

	private void updateFromModel() {
		runToolCheckBox.setSelected(model.getRunCreateIndex());
		proteinDBPanel.updateFromModel();
		massPanel.updateFromModel();
		digestPanel.updateFromModel();
		enzymePanel.updateFromModel();
		missedCleavages.setSelected(model.getAllowMissedCleavages());
	}
	
	public void itemStateChanged(final ItemEvent event) {
		logger.info("User selected index component.");
		updateFromModel();
		if (((CruxComponentButton) event.getSource()).isSelected() == true) {
			setVisible(true);
		}
		else {
			setVisible(false);
		}
	}
	
	class RunComponentChangeListener implements ItemListener {
		public void itemStateChanged(final ItemEvent event) {
			//saveButton.setEnabled(runToolCheckBox.isSelected());
		}
	}
	
	class SaveButtonListener implements ActionListener {
		public void actionPerformed(ActionEvent e) {
			boolean checked = runToolCheckBox.isSelected();
			model.setRunCreateIndex(checked);
			String proteinDatabase = proteinDBPanel.getProteinDatabase();
			if (checked && proteinDatabase == null || proteinDatabase.length() == 0) {
				model.setRunSearchForMatches(false);
				button.setSelectedToRun(false);
				JOptionPane.showMessageDialog(null, "A protein database must be specified to search for matches.");
			}
			else if (checked) {
				model.setProteinDatabase(proteinDBPanel.getProteinDatabase());
				model.setMassType(massPanel.getMassType());
				model.setDigestType(digestPanel.getDigestType());
				model.setEnzyme(enzymePanel.getEnzyme());
				checked = missedCleavages.isSelected();
				model.setAllowMissedCleavages(checked);	
				model.setRunCreateIndex(true);
				button.setSelectedToRun(true);
				setVisible(false);			
				dummyButton.setSelected(true);
			}
			else {
				model.setRunCreateIndex(false);
				button.setSelectedToRun(false);
				setVisible(false);			
				dummyButton.setSelected(true);
			}
		}
	}
	class CancelButtonListener implements ActionListener {
		public void actionPerformed(ActionEvent e) {
			setVisible(false);			
			dummyButton.setSelected(true);
		}
	}

}