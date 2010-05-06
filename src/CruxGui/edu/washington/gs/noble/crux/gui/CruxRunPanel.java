package edu.washington.gs.noble.crux.gui;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.logging.Logger;

import javax.swing.Box;
import javax.swing.BoxLayout;
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

public class CruxRunPanel extends JPanel {
	
	private static Logger logger = 
		Logger.getLogger("edu.washington.gs.noble.crux.gui");
	
	private final CruxGui cruxGui;
	private final JPanel outputPanel = new JPanel();
	private final JScrollPane scrollPane;
	private final JTextArea textArea;

	CruxRunPanel(CruxGui cruxGui) {
		this.cruxGui = cruxGui;
		setLayout(new BorderLayout());
		setBackground(Color.BLUE);
		setAlignmentX(Component.LEFT_ALIGNMENT);
		setVisible(false);
	    scrollPane = new JScrollPane(outputPanel, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
	    scrollPane.setSize(200, 400);
	    scrollPane.setBackground(Color.green);
	    textArea = new JTextArea(100, 200);
	    textArea.setBackground(Color.white);
	    //scrollPane.add(textArea); 
	    add(textArea, BorderLayout.CENTER);
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
}
