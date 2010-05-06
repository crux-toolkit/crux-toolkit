SOURCE = \
  edu/washington/gs/noble/crux/gui/CruxAnalysisModel.java \
  edu/washington/gs/noble/crux/gui/CruxAppButtonGroup.java \
  edu/washington/gs/noble/crux/gui/CruxComponentButton.java \
  edu/washington/gs/noble/crux/gui/CruxGui.java \
  edu/washington/gs/noble/crux/gui/CruxMainFrame.java \
  edu/washington/gs/noble/crux/gui/CruxParameterControl.java \
  edu/washington/gs/noble/crux/gui/CruxParameterPanel.java \
  edu/washington/gs/noble/crux/gui/CruxRunPanel.java \
  edu/washington/gs/noble/crux/gui/CruxSetupPanel.java \
  edu/washington/gs/noble/crux/gui/DigestPanel.java \
  edu/washington/gs/noble/crux/gui/EnzymePanel.java \
  edu/washington/gs/noble/crux/gui/FeatureFilePanel.java \
  edu/washington/gs/noble/crux/gui/ImagePanel.java \
  edu/washington/gs/noble/crux/gui/LengthPanel.java \
  edu/washington/gs/noble/crux/gui/MassPanel.java \
  edu/washington/gs/noble/crux/gui/MassTypePanel.java \
  edu/washington/gs/noble/crux/gui/MasterButtonPanel.java \
  edu/washington/gs/noble/crux/gui/MaxModsPanel.java \
  edu/washington/gs/noble/crux/gui/ProteinDBPanel.java \
  edu/washington/gs/noble/crux/gui/SpectraPanel.java \
  edu/washington/gs/noble/crux/gui/SpectrumChargePanel.java \
  edu/washington/gs/noble/crux/gui/SpectrumMassPanel.java

all: classes CruxGui.jar

clean:
	rm -f edu/washington/gs/noble/crux/gui/*.class
	rm -f CruxGui.jar

classes: $(SOURCE) edu/washington/gs/noble/crux/gui/schematic.png
	javac $(SOURCE)

CruxGui.jar: $(SOURCE) CruxGui.mf
	jar cmf \
    CruxGui.mf \
    CruxGui.jar \
    edu/washington/gs/noble/crux/gui/*.class \
    edu/washington/gs/noble/crux/gui/schematic.png 
