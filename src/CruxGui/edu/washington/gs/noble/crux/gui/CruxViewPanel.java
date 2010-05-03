package edu.washington.gs.noble.crux.gui;

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.util.logging.Logger;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JPanel;

/**
 * The CruxViewPanelClass presents the schematic diagram of data flow between the 
 * crux tools and the GUI elements needed to get and set the crux analysis parameters.
 *
 * @author Charles E. Grant
 *
 */
@SuppressWarnings("serial")
public class CruxViewPanel extends JPanel {
	
	private static Logger logger = 
		Logger.getLogger("edu.washington.gs.noble.crux.gui");
	
	private final CruxGui cruxGui;
	private final JPanel properties = new JPanel();
	private final CruxAppButtonGroup buttonGroup = new CruxAppButtonGroup();
	private final CruxComponentButton createIndexButton = new CruxComponentButton("Create Index");
	private final CruxComponentButton searchButton = new CruxComponentButton("Search For Matches");
	private final CruxComponentButton computeQValuesButton = new CruxComponentButton("Compute Q-values");
	private final CruxComponentButton percolatorButton = new CruxComponentButton("Percolator");
	private final CruxComponentButton qrankerButton = new CruxComponentButton("Q-ranker");
	private final CruxComponentButton dummyButton = new CruxComponentButton(""); // This button only exists to select, so that other button will be unselected.
	private final ImagePanel image;

	CruxViewPanel(CruxGui cruxGui) {
		
		this.cruxGui = cruxGui;
		setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));
		setBackground(Color.white);
		setAlignmentX(Component.LEFT_ALIGNMENT);
		
		// Group buttons
		buttonGroup.add(createIndexButton);
		buttonGroup.add(searchButton);
		buttonGroup.add(computeQValuesButton);
		buttonGroup.add(percolatorButton);
		buttonGroup.add(qrankerButton);
		buttonGroup.add(dummyButton);
		
		// Set up panel containing flow chart and buttons
		Component component2 = Box.createVerticalStrut(20);
		add(component2);
		image = new ImagePanel(createIndexButton, searchButton, computeQValuesButton, percolatorButton, qrankerButton);
		setSize(new Dimension(image.getWidth() + 25, image.getHeight() + 25));
		add(image);
		
		// Set up panel where controls for setting parameters will appear
		properties.setBackground(Color.white);
		properties.setSize(new Dimension(image.getWidth(), 4 * image.getHeight()));
		properties.setLayout(null);
		initIndexPanel();
		initSearchPanel();
		initQValuesPanel();
		initPercolatorPanel();
		initQRankerPanel();
		add(properties);
	}
	
	private void initIndexPanel() {
		ParameterPanel indexParameterPanel = new ParameterPanel(CruxAnalysisModel.CruxComponents.CREATE_INDEX, cruxGui, "Index creation paramters", createIndexButton, dummyButton);
		indexParameterPanel.setBounds(createIndexButton.getX()-25, 30, 400, 350);
		indexParameterPanel.addParameterControl(new ProteinDBPanel(cruxGui));
		indexParameterPanel.addParameterControl(new EnzymePanel(cruxGui));
		indexParameterPanel.addParameterControl(new MassPanel(cruxGui));
		indexParameterPanel.addParameterControl(new DigestPanel(cruxGui));
		properties.add(indexParameterPanel);
	}
	
	private void initSearchPanel() {
		ParameterPanel searchParameterPanel = new ParameterPanel(CruxAnalysisModel.CruxComponents.SEARCH_FOR_MATCHES, cruxGui, "Search paramters", searchButton, dummyButton);
		searchParameterPanel.setBounds(searchButton.getX()-25, 30, 400, 350);
		searchParameterPanel.addParameterControl(new ProteinDBPanel(cruxGui));
		searchParameterPanel.addParameterControl(new SpectraPanel(cruxGui));
		searchParameterPanel.addParameterControl(new EnzymePanel(cruxGui));
		searchParameterPanel.addParameterControl(new MassPanel(cruxGui));
		searchParameterPanel.addParameterControl(new DigestPanel(cruxGui));
		properties.add(searchParameterPanel);
	}
	
	private void initQValuesPanel() {
		ParameterPanel qvaluesParameterPanel = new ParameterPanel(CruxAnalysisModel.CruxComponents.COMPUTE_Q_VALUES, cruxGui, "Compute q-value paramters", computeQValuesButton, dummyButton);
		qvaluesParameterPanel.setBounds(computeQValuesButton.getX()-25, 30, 400, 350);
		properties.add(qvaluesParameterPanel);
	}
	
	private void initPercolatorPanel() {
		ParameterPanel percolatorParameterPanel = new ParameterPanel(CruxAnalysisModel.CruxComponents.PERCOLATOR, cruxGui, "Percolator paramters", percolatorButton, dummyButton);
		percolatorParameterPanel.setBounds(computeQValuesButton.getX()-25, 30, 400, 350); // Line up to the leftmost of the buttons
		properties.add(percolatorParameterPanel);
	}
	
	private void initQRankerPanel() {
		ParameterPanel qrankerParameterPanel = new ParameterPanel(CruxAnalysisModel.CruxComponents.QRANKER, cruxGui, "Q-ranker paramters", qrankerButton, dummyButton);
		qrankerParameterPanel.setBounds(computeQValuesButton.getX()-25, 30, 400, 350);
		properties.add(qrankerParameterPanel);
	}
}
