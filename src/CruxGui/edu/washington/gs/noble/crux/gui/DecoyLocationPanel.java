package edu.washington.gs.noble.crux.gui;

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.util.logging.Logger;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;


@SuppressWarnings("serial")
    class DecoyLocationPanel extends CruxParameterControl {
	private static Logger logger = 
	    Logger.getLogger("edu.washington.gs.noble.crux.gui");
	private final CruxGui cruxGui;
	private final JLabel decoyLocationLabel = new JLabel("Decoy Location:");
	private final JComboBox decoyLocationCombo = 
	    new JComboBox(CruxAnalysisModel.DecoyLocation.values());
	
	public DecoyLocationPanel(CruxGui cruxGui){
	    super();
	    this.cruxGui = cruxGui;
	    this.setAlignmentX(Component.LEFT_ALIGNMENT);
	    this.setLayout(new BoxLayout(this, BoxLayout.X_AXIS));
	    CruxAnalysisModel model = cruxGui.getAnalysisModel();
	    add(decoyLocationLabel);
	    add(Box.createRigidArea(new Dimension(10,0)));
	    decoyLocationCombo.setSelectedIndex(model.getDecoyLocationDefault().ordinal());
	    decoyLocationCombo.setBackground(Color.white);
	    decoyLocationCombo.setMaximumSize(new Dimension(230, 24));
	    add(decoyLocationCombo);
	}
	
	public void loadDefaults(){
	    CruxAnalysisModel model = cruxGui.getAnalysisModel();
	    decoyLocationCombo.setSelectedIndex(model.getDecoyLocationDefault().ordinal());
	    logger.info("DecoyLocation read from model: "+model.getDecoyLocationDefault().toString());
	}
	
	public void saveToModel(){
	    CruxAnalysisModel model = cruxGui.getAnalysisModel();
	    model.setDecoyLocation(getDecoyLocation());
	}
	
	public void updateFromModel(){
	    CruxAnalysisModel model = cruxGui.getAnalysisModel();
	    decoyLocationCombo.setSelectedIndex(model.getDecoyLocation().ordinal());
	    logger.info("DecoyLocation read from model: "+model.getDecoyLocation().toString());
	}
	
	public CruxAnalysisModel.DecoyLocation getDecoyLocation() {
	    return (CruxAnalysisModel.DecoyLocation) decoyLocationCombo.getSelectedItem();
	}

    }