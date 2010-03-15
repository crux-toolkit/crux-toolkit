package edu.washington.gs.noble.crux.gui;

import java.awt.EventQueue;
import java.util.logging.Logger;

public class CruxGui implements Runnable {
	
	CruxAnalysisModel model;

	private static Logger logger = Logger
			.getLogger("edu.washington.gs.noble.crux.gui");
	
	public static void main(String[] args) {
		EventQueue.invokeLater(new CruxGui());		
	}
		
	public void run() {
		logger.info("Starting the Crux GUI");
		model = new CruxAnalysisModel();
		CruxMainFrame frame = new CruxMainFrame(this);
		frame.setVisible(true);
	}
	
	public void shutdown() {
		logger.info("Shutting down the Crux GUI");
		System.exit(0);
	}
	
	public CruxAnalysisModel getAnalysisModel() {
		return model;
	}
	
	public void setAnalysisModel(CruxAnalysisModel model) {
		this.model = model;
	}
}


