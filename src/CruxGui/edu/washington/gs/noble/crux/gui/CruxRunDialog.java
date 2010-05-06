package edu.washington.gs.noble.crux.gui;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.logging.Logger;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;


import edu.washington.gs.noble.crux.gui.CruxAnalysisModel.CruxComponents;

/**
 * This panel is used to execute the analysis and report the output 
 * to the user.
 * 
 * @author Charles E. Grant
 *
 */

public class CruxRunDialog extends JDialog {
	
	private static Logger logger = 
		Logger.getLogger("edu.washington.gs.noble.crux.gui");
	
	private final CruxGui cruxGui;
	private final JPanel outputPanel = new JPanel();
	private final JScrollPane scrollPane;
	private final JTextArea textArea;
	private final JButton closeButton = new JButton("Close Button");

	CruxRunDialog(CruxGui cruxGui) {
		super(cruxGui.frame, "Run Analysis", JDialog.DEFAULT_MODALITY_TYPE.APPLICATION_MODAL);
		this.cruxGui = cruxGui;
		setBackground(Color.BLUE);
	    setSize(200, 400);
	    closeButton.addActionListener(new CloseButtonListener());
	    add(closeButton, BorderLayout.SOUTH);
	    textArea = new JTextArea(100, 200);
	    textArea.setBackground(Color.YELLOW);
	    scrollPane = new JScrollPane(textArea, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
	    scrollPane.setSize(200, 400);
	    scrollPane.setBackground(Color.GREEN);
	    scrollPane.add(textArea); 
	    add(scrollPane, BorderLayout.CENTER);
	    addWindowListener(new RunWindowListener());
	    runAnalysis();
	}
	
	public void runAnalysis() {
		
		CruxAnalysisModel model = cruxGui.getAnalysisModel();
		boolean continueRunAnalysis = true;
		if (model != null && model.needsSaving()) {
			continueRunAnalysis = cruxGui.promptToSave();
		}
		if (continueRunAnalysis) {
			if (model.isValidAnalysis()) {
				if (model.getPathToCrux() == null) {
					continueRunAnalysis = cruxGui.promptForCruxPath();
				}
				if (continueRunAnalysis) {
			        Process process = null;
					// Write out a parameter file to the analysis directory
					if (!model.saveModelToParameterFile()) {
						logger.info("Unable to write parameter file.");
						return;
					}

			        cruxGui.frame.repaint();
			        for (CruxComponents component: CruxComponents.values()) {
			        	repaint();
    			        process = model.runComponent(component);
    			        if (process != null) {
    			        	textArea.append("Running " + component.toString() + "\n");
    			        	InputStream errStream = process.getErrorStream();
    			        	BufferedReader errReader = new BufferedReader(new InputStreamReader(errStream));
    			        	String errString;
    			        	try {
								while ((errString = errReader.readLine()) != null){ 
									textArea.append(errString + "\n");
								}
							} catch (IOException e1) {
								// TODO Auto-generated catch block
								e1.printStackTrace();
							}
	    			        try {
								process.waitFor();
							} catch (InterruptedException e) {
								// TODO Auto-generated catch block
								e.printStackTrace();
							}
							logger.info("Execution of " + component.toString() + " completed with exit status " + process.exitValue());
							if (process.exitValue() == 0) {
								model.setRunStatus(component, CruxAnalysisModel.RunStatus.COMPLETED);
							}
							else {
								model.setRunStatus(component, CruxAnalysisModel.RunStatus.FAILED);
								break;
							}
    			        }
			        }
			        cruxGui.frame.updateFromModel(model);
			        cruxGui.frame.repaint();
				}
			}
		}
		
	}
	
	class CloseButtonListener implements ActionListener {
		
		public void actionPerformed(final ActionEvent e) {
			setVisible(false);
		}
		
	}
	
	class RunWindowListener implements WindowListener {

		/* (non-Javadoc)
		 * @see java.awt.event.WindowListener#windowActivated(java.awt.event.WindowEvent)
		 */
		@Override
		public void windowActivated(WindowEvent e) {
			// TODO Auto-generated method stub
			
		}

		/* (non-Javadoc)
		 * @see java.awt.event.WindowListener#windowClosed(java.awt.event.WindowEvent)
		 */
		@Override
		public void windowClosed(WindowEvent e) {
			// TODO Auto-generated method stub
			
		}

		/* (non-Javadoc)
		 * @see java.awt.event.WindowListener#windowClosing(java.awt.event.WindowEvent)
		 */
		@Override
		public void windowClosing(WindowEvent e) {
			// TODO Auto-generated method stub
			
		}

		/* (non-Javadoc)
		 * @see java.awt.event.WindowListener#windowDeactivated(java.awt.event.WindowEvent)
		 */
		@Override
		public void windowDeactivated(WindowEvent e) {
			// TODO Auto-generated method stub
			
		}

		/* (non-Javadoc)
		 * @see java.awt.event.WindowListener#windowDeiconified(java.awt.event.WindowEvent)
		 */
		@Override
		public void windowDeiconified(WindowEvent e) {
			// TODO Auto-generated method stub
			
		}

		/* (non-Javadoc)
		 * @see java.awt.event.WindowListener#windowIconified(java.awt.event.WindowEvent)
		 */
		@Override
		public void windowIconified(WindowEvent e) {
			// TODO Auto-generated method stub
			
		}

		/* (non-Javadoc)
		 * @see java.awt.event.WindowListener#windowOpened(java.awt.event.WindowEvent)
		 */
		@Override
		public void windowOpened(WindowEvent e) {
			runAnalysis();
			
		}
	
	}

}