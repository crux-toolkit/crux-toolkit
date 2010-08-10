package edu.washington.gs.noble.crux.gui;

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.util.logging.Logger;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JPanel;

/**
 * The CruxSetupPanelClass presents the schematic diagram of data flow between the 
 * crux tools and the GUI elements needed to get and set the crux analysis parameters.
 *
 * @author Charles E. Grant
 *
 */
@SuppressWarnings("serial")
public class CruxSetupPanel extends JPanel {
	
	private static Logger logger = 
		Logger.getLogger("edu.washington.gs.noble.crux.gui");
	
	private final CruxGui cruxGui;
	private final JPanel properties = new JPanel();
	private final CruxAppButtonGroup buttonGroup = new CruxAppButtonGroup();
	private final CruxComponentButton createIndexButton = new CruxComponentButton("Create Index", CruxAnalysisModel.CruxComponents.CREATE_INDEX);
	private final CruxComponentButton searchButton = new CruxComponentButton("Search For Matches", CruxAnalysisModel.CruxComponents.SEARCH_FOR_MATCHES);
	private final CruxComponentButton computeQValuesButton = new CruxComponentButton("Compute Q-values", CruxAnalysisModel.CruxComponents.COMPUTE_Q_VALUES);
	private final CruxComponentButton percolatorButton = new CruxComponentButton("Percolator", CruxAnalysisModel.CruxComponents.PERCOLATOR);
	private final CruxComponentButton qrankerButton = new CruxComponentButton("Q-ranker", CruxAnalysisModel.CruxComponents.QRANKER);
	private final CruxComponentButton dummyButton = new CruxComponentButton("", null); // This button only exists to select, so that other button will be unselected.
        private CruxParameterPanel indexParameterPanel, searchParameterPanel, qvaluesParameterPanel, percolatorParameterPanel,qrankerParameterPanel;
	private final ImagePanel image;

	CruxSetupPanel(CruxGui cruxGui) {
		
		this.cruxGui = cruxGui;
		setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));
		setBackground(Color.white);
		setAlignmentX(Component.LEFT_ALIGNMENT);
		createIndexButton.requestFocusInWindow();
		
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
		indexParameterPanel = new CruxParameterPanel(CruxAnalysisModel.CruxComponents.CREATE_INDEX, cruxGui, "Index creation paramters", createIndexButton, dummyButton);
		indexParameterPanel.setBounds(createIndexButton.getX()-25, 30, 400, 350);
		indexParameterPanel.addParameterControl(new ProteinDBPanel(cruxGui));
		indexParameterPanel.addParameterControl(new MissedCleavagesPanel(cruxGui));
		indexParameterPanel.addParameterControl(new EnzymePanel(cruxGui));
		indexParameterPanel.addParameterControl(new MassTypePanel(cruxGui));
		indexParameterPanel.addParameterControl(new DigestPanel(cruxGui));
		indexParameterPanel.addParameterControl(new MassPanel(cruxGui));
		indexParameterPanel.addParameterControl(new LengthPanel(cruxGui));
		indexParameterPanel.addParameterControl(new VerbosityPanel(cruxGui));
		indexParameterPanel.addParameterControl(new CustomEnzymeCleavagePanel(cruxGui));
		
		properties.add(indexParameterPanel);
	}
	
	private void initSearchPanel() {
		searchParameterPanel = new CruxParameterPanel(CruxAnalysisModel.CruxComponents.SEARCH_FOR_MATCHES, cruxGui, "Search paramters", searchButton, dummyButton);
		searchParameterPanel.setBounds(searchButton.getX()-25, 30, 400, 350);
		searchParameterPanel.addParameterControl(new DecoyLocationPanel(cruxGui));
		searchParameterPanel.addParameterControl(new SpectraPanel(cruxGui));
		searchParameterPanel.addParameterControl(new NumDecoysPerTargetPanel(cruxGui));
		searchParameterPanel.addParameterControl(new SpectrumMassPanel(cruxGui));
		searchParameterPanel.addParameterControl(new SpectrumChargePanel(cruxGui));
		searchParameterPanel.addParameterControl(new MaxModsPanel(cruxGui));
		searchParameterPanel.addParameterControl(new VerbosityPanel(cruxGui));
		properties.add(searchParameterPanel);
	}
	
	private void initQValuesPanel() {
		qvaluesParameterPanel = new CruxParameterPanel(CruxAnalysisModel.CruxComponents.COMPUTE_Q_VALUES, cruxGui, "Compute q-value paramters", computeQValuesButton, dummyButton);
		qvaluesParameterPanel.addParameterControl(new VerbosityPanel(cruxGui)); 
		qvaluesParameterPanel.setBounds(computeQValuesButton.getX()-60, 30, 400, 350);
		properties.add(qvaluesParameterPanel);
	}
	
	private void initPercolatorPanel() {
		percolatorParameterPanel = new CruxParameterPanel(CruxAnalysisModel.CruxComponents.PERCOLATOR, cruxGui, "Percolator paramters", percolatorButton, dummyButton);
		percolatorParameterPanel.setBounds(computeQValuesButton.getX()-60, 30, 400, 350); // Line up to the leftmost of the buttons
		percolatorParameterPanel.addParameterControl(new FeatureFilePanel(cruxGui));
		percolatorParameterPanel.addParameterControl(new VerbosityPanel(cruxGui));
		properties.add(percolatorParameterPanel);
	}
	
	private void initQRankerPanel() {
		qrankerParameterPanel = new CruxParameterPanel(CruxAnalysisModel.CruxComponents.QRANKER, cruxGui, "Q-ranker paramters", qrankerButton, dummyButton);
		qrankerParameterPanel.setBounds(computeQValuesButton.getX()-60, 30, 400, 350);
		qrankerParameterPanel.addParameterControl(new FeatureFilePanel(cruxGui));
		qrankerParameterPanel.addParameterControl(new VerbosityPanel(cruxGui));
		properties.add(qrankerParameterPanel);
	}

    public void setInitialFocus(){
    	createIndexButton.requestFocusInWindow();
    }

	
    public void updateFromModel(CruxAnalysisModel model) {
    	createIndexButton.updateFromModel(model);
    	searchButton.updateFromModel(model);
    	computeQValuesButton.updateFromModel(model);
    	percolatorButton.updateFromModel(model);
    	qrankerButton.updateFromModel(model);

	indexParameterPanel.updateFromModel();
	searchParameterPanel.updateFromModel();
	qvaluesParameterPanel.updateFromModel();
	percolatorParameterPanel.updateFromModel();
	qrankerParameterPanel.updateFromModel();

	
    }
}
