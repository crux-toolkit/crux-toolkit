package edu.washington.gs.noble.crux.gui;

import javax.swing.JPanel;

/**
 * This is an abstract class used to extend the properties of the Swing 
 * JPanel class and to require the implementation of the functions:
 * 	  updateFromModel()
 * 	  saveToModel()
 * 
 * This allows to create several concrete classes representing GUI elements 
 * for settting parameters for a crux analysis. Each such class must be able to
 * load and save its state from the current CruxAnalysisModel instance.
 * 
 * @author Charles E. Grant
 * 
 */
@SuppressWarnings("serial")
public abstract class CruxParameterControl extends JPanel {
	public void loadDefaults(){};
	public void updateFromModel(){};
	public void saveToModel(){};
}
