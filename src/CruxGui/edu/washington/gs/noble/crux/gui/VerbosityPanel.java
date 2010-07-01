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
    class VerbosityPanel extends CruxAdvancedParameterControl {
	private static Logger logger = 
	    Logger.getLogger("edu.washington.gs.noble.crux.gui");
	private final CruxGui cruxGui;
	private final JLabel verbosityLabel = new JLabel("verbosity:");
	private final JComboBox verbosityCombo = 
	    new JComboBox(CruxAnalysisModel.Verbosity.values());
	
	public VerbosityPanel(CruxGui cruxGui){
	    super();
	    this.cruxGui = cruxGui;
	    this.setAlignmentX(Component.LEFT_ALIGNMENT);
	    this.setLayout(new BoxLayout(this, BoxLayout.X_AXIS));
	    add(verbosityLabel);
	    add(Box.createRigidArea(new Dimension(10,0)));
	    verbosityCombo.setSelectedIndex(0);
	    verbosityCombo.setBackground(Color.white);
	    verbosityCombo.setMaximumSize(new Dimension(250, 24));
	    add(verbosityCombo);
	}
	
	public void loadDefaults(){
	    CruxAnalysisModel model = cruxGui.getAnalysisModel();
	    verbosityCombo.setSelectedIndex(model.getVerbosityDefault().ordinal());
	    logger.info("Verbosity read from model: "+model.getVerbosityDefault().toString());
	}
	
	public void saveToModel(){
	    CruxAnalysisModel model = cruxGui.getAnalysisModel();
	    model.setVerbosity(getVerbosity());
	}
	
	public void updateFromModel(){
	    CruxAnalysisModel model = cruxGui.getAnalysisModel();
	    verbosityCombo.setSelectedIndex(model.getVerbosity().ordinal());
	    logger.info("Verbosity read from model: "+model.getVerbosity().toString());
	}
	
	public CruxAnalysisModel.Verbosity getVerbosity() {
	    return (CruxAnalysisModel.Verbosity) verbosityCombo.getSelectedItem();
	}

    }