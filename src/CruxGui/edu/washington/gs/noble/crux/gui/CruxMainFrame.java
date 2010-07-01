package edu.washington.gs.noble.crux.gui;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.util.logging.Logger;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;

/**
 * This class specifies the main window (JFrame) for the CruxGui.
 * There are three sub-components: a CruxSetupPanel, a CruxRunPanel, and a MasterButton panel.
 * 
 * The CruxSetupPanel contains the a schmatic diagram of the flow of the the
 * crux tools, overlayed with buttons, and panels containing controls for setting
 * the parameters of a crux analysis.
 * 
 * The CruxRunPanel displays the standard error output from the analysis
 * 
 * The MasterButton panel contains several buttons for naming, loading, and running
 * analysis, and for exiting the application.
 * 
 * @author Charles E. Grant
 * 
 */

@SuppressWarnings("serial")
public class CruxMainFrame extends JFrame {

	private static Logger logger = Logger
		.getLogger("edu.washington.gs.noble.crux.gui");

	private final CruxGui cruxGui;
	public final CruxSetupPanel cruxSetupPanel;
	public final CruxRunPanel cruxRunPanel;
	private final MasterButtonPanel buttonPanel;
	
	public  CruxMainFrame(final CruxGui cruxGui) {
		this.cruxGui = cruxGui;
		cruxSetupPanel = new CruxSetupPanel(cruxGui);
		cruxRunPanel = new CruxRunPanel(cruxGui);
		buttonPanel = new MasterButtonPanel(cruxGui);
		setResizable(false);
		setUndecorated(false);
		setTitle("Crux");
		setBackground(Color.WHITE);
		setLayout(new BorderLayout());	
		addWindowListener(new CruxWindowAdapter());
		add(buttonPanel, BorderLayout.SOUTH);
		add(cruxSetupPanel, BorderLayout.CENTER);
		setSize(cruxSetupPanel.getWidth(), 4 * cruxSetupPanel.getHeight());
		setLocationRelativeTo(null);
	}
	
	public void updateFromModel(CruxAnalysisModel model) {
		cruxSetupPanel.updateFromModel(model);
		String name = model.getName();
		if (name != null) {
		    setTitle("Crux - " + name);
		}
		else {
		    setTitle("Crux");
		}
		repaint();
	}
	
	private class CruxWindowAdapter extends WindowAdapter {
		public void windowClosing(WindowEvent we) {
			cruxGui.shutdown();
		}
		
	}
    
    public void setInitialFocus(){
	//buttonPanel.setInitialFocus();
	cruxSetupPanel.setInitialFocus();
    }

}