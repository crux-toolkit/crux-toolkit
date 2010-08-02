package edu.washington.gs.noble.crux.gui;

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.Vector;
import java.util.logging.Logger;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.border.EtchedBorder;

@SuppressWarnings("serial")
class CruxParameterPanel extends JPanel implements ItemListener {
	
	private static Logger logger = 
		Logger.getLogger("edu.washington.gs.noble.crux.gui");
	
	private final CruxGui cruxGui;
	private final Vector<CruxParameterControl> parameterControls = new Vector<CruxParameterControl>();
	private final JPanel panel = new JPanel();
	private final JCheckBox runToolCheckBox = new JCheckBox("Run this tool");
	private final JCheckBox showAdvancedParameters = new JCheckBox("Show advanced parameters");
	private final JButton saveButton = new JButton("Save");
	private final JButton cancelButton = new JButton("Cancel");
	private final JButton loadDefaultsButton = new JButton("Load Defaults");
	private final CruxAnalysisModel.CruxComponents component;
	private final JScrollPane scrollPane;
	private final CruxComponentButton button;
	private final CruxComponentButton dummyButton;
	
	public CruxParameterPanel(CruxAnalysisModel.CruxComponents component, CruxGui cruxGui, final String title, final CruxComponentButton button, final CruxComponentButton dummyButton) {
		super();
		this.cruxGui = cruxGui;
		this.component = component;
		this.button = button;
		this.button.addItemListener(this);
		this.dummyButton = dummyButton;
		setBorder(
			BorderFactory.createTitledBorder(
				BorderFactory.createEtchedBorder(EtchedBorder.RAISED), 
				title
			)
		);
		
		// Top controls
		setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
		runToolCheckBox.addItemListener(new RunComponentChangeListener());
		showAdvancedParameters.addActionListener(new showAdvancedParametersListener());

		add(runToolCheckBox);
		add(Box.createRigidArea(new Dimension(0,12)));
		add(showAdvancedParameters);
		add(Box.createRigidArea(new Dimension(0,12)));
		
		// Scroll pane controls
	    scrollPane = new JScrollPane(panel, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
	    add(scrollPane);
		panel.setLayout(new BoxLayout(panel, BoxLayout.Y_AXIS));
		panel.add(Box.createRigidArea(new Dimension(0,12)));
		add(Box.createRigidArea(new Dimension(0,12)));
	    
	    // Bottom button panel controls
		JPanel buttonPanel = new JPanel();
		buttonPanel.setLayout(new BoxLayout(buttonPanel, BoxLayout.LINE_AXIS));
		buttonPanel.setAlignmentX(Component.LEFT_ALIGNMENT);
		saveButton.addActionListener(new SaveButtonListener());
		buttonPanel.add(Box.createRigidArea(new Dimension(36,0)));
		buttonPanel.add(saveButton);
		cancelButton.addActionListener(new CancelButtonListener());
		buttonPanel.add(Box.createRigidArea(new Dimension(12,0)));
		buttonPanel.add(cancelButton);
		buttonPanel.add(Box.createRigidArea(new Dimension(12,0)));
		loadDefaultsButton.addActionListener(new LoadDefaultsButtonListener());
		buttonPanel.add(loadDefaultsButton);
		add(buttonPanel);
		
		updateFromModel();
		setVisible(false);
	}
	
	public void addParameterControl(CruxParameterControl control) {
		parameterControls.add(control);
		panel.add(control);
		panel.add(Box.createRigidArea(new Dimension(0, 20)));
	}
	
	public void addParameterControls(Vector<CruxParameterControl> controls) {
		for (CruxParameterControl control: controls) {
			parameterControls.add(control);
			panel.add(control);
			panel.add(Box.createRigidArea(new Dimension(0, 20)));
		}
	}
    
    public void updateFromModel() {
		CruxAnalysisModel model = cruxGui.getAnalysisModel();
		runToolCheckBox.setSelected(model.getRunComponent(component));
		showAdvancedParameters.setSelected(model.getShowAdvancedParameters(component));
		for (CruxParameterControl component: parameterControls) {
			component.updateFromModel();
		}
	}
	
	public void itemStateChanged(final ItemEvent event) {
		logger.info("User selected index component.");
		CruxAnalysisModel model = cruxGui.getAnalysisModel();
		saveState(model);
		updateFromModel();
		showAdvancedParameters();
		if (((CruxComponentButton) event.getSource()).isSelected() == true) {
			setVisible(true);
		}
		else {
			setVisible(false);
		}
		
	}

    public void showAdvancedParameters(){
	boolean show = showAdvancedParameters.isSelected();
	for (CruxParameterControl component: parameterControls){
	    if (component instanceof CruxAdvancedParameterControl){
		component.setVisible(show);
	    }
	}
    }

    public void saveState(CruxAnalysisModel model){
	boolean checked = runToolCheckBox.isSelected();
	model.setRunComponent(component, checked);
	if (checked) {
	    button.setSelectedToRun(true);
	    model.setRunComponent(component, true);
	}
	else {
	    button.setSelectedToRun(false);
	    model.setRunComponent(component, false);
	}
	for (CruxParameterControl control: parameterControls) {
	    control.saveToModel();
	}
	model.setRunStatus(component, CruxAnalysisModel.RunStatus.NOT_RUN);
	model.setNeedsSaving(true);
    }

    
	
	class RunComponentChangeListener implements ItemListener {
		public void itemStateChanged(final ItemEvent event) {
		}
	}
	
class LoadDefaultsButtonListener implements ActionListener {
		public void actionPerformed(ActionEvent e) {
			CruxAnalysisModel model = cruxGui.getAnalysisModel();
			for (CruxParameterControl component: parameterControls) {
				component.loadDefaults();
			}
		}
	}
	
	class SaveButtonListener implements ActionListener {
		public void actionPerformed(ActionEvent e) {
		        CruxAnalysisModel model = cruxGui.getAnalysisModel();
		        saveState(model);
			setVisible(false);			
			dummyButton.setSelected(true);
			cruxGui.frame.updateFromModel(model);
		}
	}
	
	class CancelButtonListener implements ActionListener {
	    public void actionPerformed(ActionEvent e) {
		   CruxAnalysisModel model = cruxGui.getAnalysisModel();
		   updateFromModel();
		    cruxGui.frame.updateFromModel(model);
			setVisible(false);			
			dummyButton.setSelected(true);
		}
	}

    class showAdvancedParametersListener implements ActionListener {
	public void actionPerformed(ActionEvent e){
	    showAdvancedParameters();
	}
    }
}