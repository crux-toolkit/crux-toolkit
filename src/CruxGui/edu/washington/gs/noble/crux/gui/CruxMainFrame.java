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

@SuppressWarnings("serial")
public class CruxMainFrame extends JFrame {

	private static Logger logger = Logger
		.getLogger("edu.washington.gs.noble.crux.gui");

	CruxAnalysisModel model;
	public CruxViewPanel cruxGuiPanel;
	MasterButtonPanel buttonPanel;
	
	public  CruxMainFrame(final CruxGui gui) {
		setResizable(false);
		setUndecorated(false);
		setTitle("Crux");
		setBackground(Color.white);
		setLayout(new BorderLayout());	
		addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent we) {
				gui.shutdown();
			}
		});
		model = gui.getAnalysisModel();
		cruxGuiPanel = new CruxViewPanel(model);
		buttonPanel = new MasterButtonPanel(gui);
		add(buttonPanel, BorderLayout.SOUTH);
		add(cruxGuiPanel, BorderLayout.CENTER);
		setSize(cruxGuiPanel.getWidth(), 4 * cruxGuiPanel.getHeight());
		setLocationRelativeTo(null);
	}

}