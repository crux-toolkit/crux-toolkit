package edu.washington.gs.noble.crux.gui;

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.util.logging.Logger;

import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JToggleButton;
import javax.swing.JPanel;

@SuppressWarnings("serial")
public class CruxViewPanel extends JPanel {
	
	private static Logger logger = 
		Logger.getLogger("edu.washington.gs.noble.crux.gui");
	
	ImagePanel image;
	final JPanel properties = new JPanel();
	final CruxAppButtonGroup buttonGroup = new CruxAppButtonGroup();
	final CruxComponentButton createIndex = new CruxComponentButton("create index");
	final CruxComponentButton search = new CruxComponentButton("search for matches");
	final CruxComponentButton computeQvalues = new CruxComponentButton("compute-q-values");
	final CruxComponentButton percolator = new CruxComponentButton("percolator");
	final CruxComponentButton qranker = new CruxComponentButton("q-ranker");
	final JToggleButton dummy = new JToggleButton(); // Only exists to select, so other button are unselected.

	CruxViewPanel(CruxAnalysisModel model) {
		
		setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));
		setBackground(Color.white);
		setAlignmentX(Component.LEFT_ALIGNMENT);
		
		// Group buttons
		buttonGroup.add(createIndex);
		buttonGroup.add(search);
		buttonGroup.add(computeQvalues);
		buttonGroup.add(percolator);
		buttonGroup.add(qranker);
		buttonGroup.add(dummy);
		
		// Set up panel with flow chart and buttons
		image = new ImagePanel(createIndex, search, computeQvalues, percolator, qranker);
		setSize(new Dimension(image.getWidth(), image.getHeight()));
		add(image);
		
		// Setup panel where controls for setting parameters will appear
		properties.setBackground(Color.white);
		properties.setSize(new Dimension(image.getWidth(), 4 * image.getHeight()));
		properties.setLayout(null);
		IndexPanel indexPanel = new IndexPanel(model, createIndex, dummy);
		indexPanel.setBounds(createIndex.getX(), 0, indexPanel.getMaximumSize().width, indexPanel.getMaximumSize().height);
		properties.add(indexPanel);
		SearchPanel searchPanel = new SearchPanel(model, search, dummy);
		searchPanel.setBounds(search.getX(), 0, searchPanel.getMaximumSize().width, searchPanel.getMaximumSize().height);
		properties.add(indexPanel);
		properties.add(searchPanel);
		add(properties);
		
		// Connect applications buttons to option panels
		createIndex.addItemListener(indexPanel);
		search.addItemListener(searchPanel);
		
	}
}
