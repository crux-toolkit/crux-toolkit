#include "LocalizeModification.h"
#include "io/MatchCollectionParser.h"
#include "tide/max_mz.h"
#include "tide/modifications.h"
#include "util/FileUtils.h"
#include "util/Params.h"
#include "io/SpectrumCollectionFactory.h"
#include "io/SpectrumRecordSpectrumCollection.h"
#include "io/SpectrumRecordWriter.h"
#include "util/StringUtils.h"
#include "TideSearchApplication.h"

using namespace std;

LocalizeModificationApplication::LocalizeModificationApplication() {
  for (int i = 0; i <= 100; i++) {
    progress_.insert(i);
  }
}

LocalizeModificationApplication::~LocalizeModificationApplication() {
}

/* The main localize-modification function.
 * Iterates over each match in a PSM results file twice:
 *   1. Look at spectrum files for each match
 *     - Check that the file exists
 *     - Load spectra into memory as a SpectrumCollection
 *   2. Search modified peptides against spectrum
 *     - For each PSM, generate a modified version of the peptide for each residue
 *       where the modification is (spectrum neutral mass - peptide mass)
 *     - Score each of these modified peptides against the spectrum
 *     - Report to output file
 */
int LocalizeModificationApplication::main(int argc, char** argv) {
  string inputFile = Params::GetString("input PSM file");
  carp(CARP_INFO, "Parsing %s...", inputFile.c_str());
  MatchCollectionParser parser;
  MatchCollection* matches = parser.create(Params::GetString("input PSM file"), "");

  map<string, string> spectrumFiles;

  uint64_t numSteps = 0;
  bool hasTargets = false;
  bool hasDecoys = false;

  MatchIterator* matchIter = new MatchIterator(matches);
  while (matchIter->hasNext()) {
    Crux::Match* match = matchIter->next();
    numSteps += match->getPeptide()->getLength() + 1;
    if (!hasTargets && !match->getNullPeptide()) {
      hasTargets = true;
    }
    if (!hasDecoys && match->getNullPeptide()) {
      hasDecoys = true;
    }
    string matchPath = match->getFilePath();
    if (spectrumFiles.find(matchPath) != spectrumFiles.end()) {
      continue;
    }
    string spectrumFile = matchPath;
    if (!FileUtils::Exists(spectrumFile)) {
      string basename = FileUtils::BaseName(spectrumFile);
      if (!FileUtils::Exists(basename)) {
        carp(CARP_FATAL, "Spectrum file '%s' could not be found", spectrumFile.c_str());
      }
      spectrumFile = basename;
    }
    spectrumFiles[matchPath] = spectrumFile;
  }
  delete matchIter;

  string outpath;
  if (hasTargets && hasDecoys) {
    outpath = make_file_path(getName() + ".txt");
  } else if (hasTargets) {
    outpath = make_file_path(getName() + ".target.txt");
  } else {
    outpath = make_file_path(getName() + ".decoy.txt");
  }

  map<string, Crux::SpectrumCollection*> spectrumCollections;
  for (map<string, string>::const_iterator i = spectrumFiles.begin(); i != spectrumFiles.end(); i++) {
    carp(CARP_INFO, "Parsing spectrum file %s", i->second.c_str());
    spectrumCollections[i->first] = SpectrumCollectionFactory::create(i->second);
    spectrumCollections[i->first]->parse();
  }

  carp(CARP_INFO, "Scoring modified peptides (results will be written to %s)...",
       FileUtils::BaseName(outpath).c_str());
  MatchFileWriter writer(outpath.c_str());
  writer.addColumnNames(this, hasTargets && hasDecoys);
  writer.writeHeader();

  double binWidth = Params::GetDouble("mz-bin-width");
  double binOffset = Params::GetDouble("mz-bin-offset");
  TheoreticalPeakSetBIons tps(200);
  tps.binWidth_ = binWidth;
  tps.binOffset_ = binOffset;

  int topMatch = Params::GetInt("top-match");
  uint64_t curStep = 0;
  matchIter = new MatchIterator(matches);
  while (matchIter->hasNext()) {
    Crux::Match* match = matchIter->next();
    Crux::Spectrum* cruxSpectrum = match->getSpectrum();
    int scan = match->getSpectrum()->getFirstScan();
    string spectrumFile = match->getFilePath();
    Crux::SpectrumCollection* collection = spectrumCollections[spectrumFile];
    if ((cruxSpectrum = collection->getSpectrum(scan)) == NULL) {
      carp(CARP_FATAL, "Spectrum %d not found in %s", scan, spectrumFile.c_str());
    } else if (cruxSpectrum->getNumPeaks() == 0) {
      delete cruxSpectrum;
      carp(CARP_WARNING, "Spectrum %d had 0 peaks, skipping", scan);
      continue;
    }
    cruxSpectrum->sortPeaks(_PEAK_LOCATION);
    double precursorMz = cruxSpectrum->getPrecursorMz();
    int charge = match->getCharge();
    Spectrum spectrum(scan, precursorMz);
    spectrum.AddChargeState(charge);
    spectrum.ReservePeaks(cruxSpectrum->getNumPeaks());
    for (PeakIterator i = cruxSpectrum->begin(); i != cruxSpectrum->end(); i++) {
      spectrum.AddPeak((*i)->getLocation(), (*i)->getIntensity());
    }
    delete cruxSpectrum;

    // Create proteins/peptides
    VariableModTable* modTable = getModTable(match);
    MassConstants::Init(modTable->ParsedModTable(),
                        modTable->ParsedNtpepModTable(), modTable->ParsedCtpepModTable(),
                        binWidth, binOffset);
    Crux::Peptide* cruxPeptide = match->getPeptide();
    vector<const pb::Protein*> proteins = createPbProteins(cruxPeptide);
    vector<pb::AuxLocation> auxLocs;
    vector<pb::Peptide> peptides = createPbPeptides(match, modTable, &auxLocs);

    // Score each peptide
    carp(CARP_DETAILED_INFO, "Scoring modified forms of %s against spectrum %d",
         cruxPeptide->getModifiedSequenceWithMasses().c_str(), scan);
    Results results(modTable);
    double neutralMass = match->getNeutralMass();
    int maxPrecursorMass = MassConstants::mass2bin(neutralMass + MAX_XCORR_OFFSET + 30) + 50;
    vector<double> evidence = spectrum.CreateEvidenceVector(binWidth, binOffset, charge,
      (MassConstants::mass2bin(cruxPeptide->calcModifiedMass()) - 0.5 + binOffset) * binWidth, maxPrecursorMass);
    for (vector<pb::Peptide>::const_iterator i = peptides.begin(); i != peptides.end(); i++) {
      Peptide peptide(*i, proteins);
      tps.Clear();
      peptide.ComputeBTheoreticalPeaks(&tps);

      if (i == peptides.begin() + 1) {
        // After we've scored the unmodified peptide, create new evidence vector for scoring the modified peptides
        evidence = spectrum.CreateEvidenceVector(binWidth, binOffset, charge,
          (MassConstants::mass2bin(neutralMass) - 0.5 + binOffset) * binWidth, maxPrecursorMass);
      }

      double xcorr = 0;
      for (vector<unsigned int>::const_iterator j = tps.unordered_peak_list_.begin();
          j != tps.unordered_peak_list_.end();
          j++) {
        xcorr += evidence[*j];
      }
      results.Add(cruxPeptide, &peptide, xcorr / 10000);
    }
    delete modTable;
    for (vector<const pb::Protein*>::const_iterator i = proteins.begin(); i != proteins.end(); i++) {
      delete *i;
    }

    // Write to output file
    results.Sort();
    for (size_t i = 0; i < topMatch && i < results.Size(); i++) {
      Crux::Peptide& peptide = *(results.Peptide(i));
      char* flanking = peptide.getFlankingAAs();
      string flankingStr(flanking);
      free(flanking);
      writer.setColumnCurrentRow(FILE_COL,                  match->getFilePath());
      writer.setColumnCurrentRow(SCAN_COL,                  scan);
      writer.setColumnCurrentRow(CHARGE_COL,                charge);
      writer.setColumnCurrentRow(SPECTRUM_PRECURSOR_MZ_COL, precursorMz);
      writer.setColumnCurrentRow(SPECTRUM_NEUTRAL_MASS_COL, neutralMass);
      writer.setColumnCurrentRow(PEPTIDE_MASS_COL,          peptide.calcModifiedMass());
      writer.setColumnCurrentRow(XCORR_SCORE_COL,           results.XCorr(i));
      writer.setColumnCurrentRow(SEQUENCE_COL,              peptide.getModifiedSequenceWithMasses());
      writer.setColumnCurrentRow(MODIFICATIONS_COL,         peptide.getModsString());
      writer.setColumnCurrentRow(PROTEIN_ID_COL,            peptide.getProteinIdsLocations());
      writer.setColumnCurrentRow(FLANKING_AA_COL,           flankingStr);
      writer.setColumnCurrentRow(TARGET_DECOY_COL,          match->isDecoy() ? "decoy" : "target");
      writer.writeRow();
    }

    curStep += cruxPeptide->getLength() + 1;
    reportProgress(curStep, numSteps);
  }
  delete matchIter;
  delete matches;

  for (map<string, Crux::SpectrumCollection*>::const_iterator i = spectrumCollections.begin();
       i != spectrumCollections.end();
       i++) {
    delete i->second;
  }

  return 0;
}

string LocalizeModificationApplication::getName() const { return "localize-modification"; }

string LocalizeModificationApplication::getDescription() const {
  return
    "This command finds, for each peptide-spectrum match (PSM) in a given set, "
    "the most likely location along the peptide for a post-translational "
    "modification (PTM). The mass of the PTM is inferred from the difference "
    "between the spectrum neutral mass and the peptide mass.";
}

vector<string> LocalizeModificationApplication::getArgs() const {
  string arr[] = { "input PSM file" };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

vector<string> LocalizeModificationApplication::getOptions() const {
  string arr[] = {
    "min-mod-mass",
    "mod-precision",
    "top-match",
    "output-dir",
    "overwrite",
    "parameter-file",
    "verbosity"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

vector< pair<string, string> > LocalizeModificationApplication::getOutputs() const {
  vector< pair<string, string> > outputs;
  outputs.push_back(make_pair(getName() + "[.target|.decoy].txt",
    "a tab-delimited text file containing the target PSMs. See <a href=\""
    "../file-formats/txt-format.html\">txt file format</a> for a list of the fields. "
    "The filename depends on whether the input file contained targets, decoys, or both."));
  outputs.push_back(make_pair(getName() + ".params.txt",
    "a file containing the name and value of all parameters/options for the "
    "current operation. Not all parameters in the file may have been used in "
    "the operation. The resulting file can be used with the --parameter-file "
    "option for other Crux programs."));
  outputs.push_back(make_pair(getName() + ".log.txt",
    "a log file containing a copy of all messages that were printed to the "
    "screen during execution."));
  return outputs;
}

bool LocalizeModificationApplication::needsOutputDirectory() const { return true; }

COMMAND_T LocalizeModificationApplication::getCommand() const { return LOCALIZE_MODIFICATION_COMMAND; }

bool LocalizeModificationApplication::hidden() const { return false; }

void LocalizeModificationApplication::reportProgress(uint64_t curStep, uint64_t totalSteps) {
  int percent = (int)((double)curStep / (double)totalSteps * 100.0);
  set<int>::iterator i = progress_.find(percent);
  if (i != progress_.end()) {
    progress_.erase(i);
    carp(CARP_INFO, "%d%% complete", percent);
  }
}

vector<const pb::Protein*> LocalizeModificationApplication::createPbProteins(
  Crux::Peptide* peptide
) const {
  vector<const pb::Protein*> proteins;
  int proteinId = -1;
  for (PeptideSrcIterator i = peptide->getPeptideSrcBegin();
       i != peptide->getPeptideSrcEnd();
       i++) {
    pb::Protein* protein = new pb::Protein();
    protein->set_id(++proteinId);
    Crux::Protein* cruxProtein = (*i)->getParentProtein();
    protein->set_name(cruxProtein->getId());
    int start = (*i)->getStartIdxOriginal();
    int length = (int)peptide->getLength();
    string residues(start + length, 'X');
    char* seq = peptide->getSequence();
    string sequence(seq);
    free(seq);
    residues.replace(start, length, sequence);
    if (cruxProtein->isPostProcess()) {
      if (start > 0) {
        char flankN = ((PostProcessProtein*)cruxProtein)->getNTermFlankingAA();
        if (flankN != '-') {
          residues[start-1] = flankN;
        }
      }
      char flankC = ((PostProcessProtein*)cruxProtein)->getCTermFlankingAA();
      if (flankC != '-') {
        residues.push_back(flankC);
      }
    }
    protein->set_residues(residues);
    proteins.push_back(protein);
  }
  return proteins;
}

vector<pb::Peptide> LocalizeModificationApplication::createPbPeptides(
  Crux::Match* match,
  VariableModTable* modTable,
  vector<pb::AuxLocation>* outAuxLocs
) const {
  if (outAuxLocs) {
    outAuxLocs->clear();
  }

  Crux::Peptide* cruxPeptide = match->getPeptide();

  vector<pb::Peptide> peptides;
  pb::Peptide peptide;
  int peptideId = -1;
  peptide.set_id(++peptideId);
  peptide.set_length(cruxPeptide->getLength());
  peptide.set_decoy_index(cruxPeptide->isDecoy() ? 0 : -1); // TODO need real decoy index if we want to write it to output

  char* seq = cruxPeptide->getSequence();
  string sequence(seq);
  free(seq);
  
  int proteinId = -1;
  pb::AuxLocation locations;
  for (PeptideSrcIterator i = cruxPeptide->getPeptideSrcBegin(); i != cruxPeptide->getPeptideSrcEnd(); i++) {
    int start = (*i)->getStartIdxOriginal();
    if (++proteinId == 0) {
      // first src
      peptide.mutable_first_location()->set_protein_id(proteinId);
      peptide.mutable_first_location()->set_pos(start);
    } else if (outAuxLocs) {
      // aux location
      pb::Location* location = locations.add_location();
      location->set_protein_id(proteinId);
      location->set_pos(start);
    }
  }
  if (outAuxLocs && locations.location_size() > 0) {
    peptide.set_aux_locations_index(0);
    outAuxLocs->push_back(locations);
  }

  // add existing modifications
  double pepMass = cruxPeptide->calcModifiedMass();
  vector<Crux::Modification> existingMods = cruxPeptide->getVarMods();
  set<unsigned char> existingIndices;
  for (vector<Crux::Modification>::const_iterator i = existingMods.begin();
       i != existingMods.end();
       i++) {
    int varModIdx = -1;
    char aa = sequence[i->Index()];
    MODS_SPEC_TYPE_T modType = modTypeToTide(i->Position());
    int numPoss = modTable->NumPoss(aa, modType);
    for (int j = 0; j < numPoss; j++) {
      int modIdx = modTable->PossDeltIx(aa, j, modType);
      // TODO Is there a better way to find the right modification?
      double modDelta = modTable->PossDelta(modIdx);
      if (abs(i->DeltaMass() - modDelta) < 0.1) {
        pepMass += modDelta;
        varModIdx = modIdx;
        break;
      }
    }
    if (varModIdx == -1) {
      carp(CARP_FATAL, "Couldn't find mod %.1f in VariableModTable", i->DeltaMass());
    }
    peptide.add_modifications(modTable->EncodeMod((int)i->Index(), varModIdx));
    existingIndices.insert(i->Index());
  }
  peptide.set_mass(pepMass);

  peptides.push_back(peptide); // Add unmodified peptide

  double modMass = calcModMass(match);
  carp(CARP_DETAILED_INFO, "Implied mod mass is %f", modMass);

  if (fabs(modMass) < Params::GetDouble("min-mod-mass")) {
    return peptides; // Implied modification mass too low, don't create modified peptides
  }

  peptide.set_mass(pepMass + modMass);

  int openModIdx = -1;
  int numPoss = modTable->NumPoss('A');
  for (int i = 0; i < numPoss; i++) {
    int modIdx = modTable->PossDeltIx('A', i);
    double modDelta = modTable->PossDelta(modIdx);
    // TODO Is there a better way to find the right modification?
    if (abs(modMass - modDelta) < 0.1) {
      openModIdx = modIdx;
      break;
    }
  }
  if (openModIdx == -1) {
    carp(CARP_FATAL, "Couldn't find diff mod in VariableModTable");
  }

  int* pbMod = peptide.mutable_modifications()->Add();
  for (unsigned char i = 0; i < cruxPeptide->getLength(); i++) {
    if (existingIndices.find(i) != existingIndices.end()) {
      continue; // don't modify same AA twice
    }
    peptide.set_id(++peptideId);
    *pbMod = modTable->EncodeMod((int)i, openModIdx);
    peptides.push_back(peptide);
  }

  return peptides;
}

MODS_SPEC_TYPE_T LocalizeModificationApplication::modTypeToTide(ModPosition position) {
  switch (position) {
  default:        return MOD_SPEC;
  case PEPTIDE_N: return NTPEP;
  case PEPTIDE_C: return CTPEP;
  case PROTEIN_N: return NTPRO;
  case PROTEIN_C: return CTPRO;
  }
}

double LocalizeModificationApplication::calcModMass(Crux::Match* match) {
  return match->getNeutralMass() - match->getPeptide()->calcModifiedMass();
}

VariableModTable* LocalizeModificationApplication::getModTable(
  Crux::Match* match
) const {
  const string allAminoAcids = "ACDEFGHIKLMNPQRSTVWY";
  const int modSpecPrecision = Params::GetInt("mass-precision");

  Crux::Peptide* peptide = match->getPeptide();
  vector<Crux::Modification> existingMods = peptide->getMods();
  char* seq = peptide->getSequence();
  string sequence(seq);
  free(seq);
  // maps contain: delta mass -> amino acids
  map< string, set<char> > modsStaticAny, modsVarAny;
  set<string> modsStaticPepN, modsVarPepN, modsStaticPepC, modsVarPepC;
  for (vector<Crux::Modification>::const_iterator i = existingMods.begin();
       i != existingMods.end();
       i++) {
    string massStr = i->DeltaMass() >= 0
      ? "+" + StringUtils::ToString(i->DeltaMass(), modSpecPrecision)
      : StringUtils::ToString(i->DeltaMass(), modSpecPrecision);
    switch (i->Position()) {
    default:
      if (i->Static()) {
        modsStaticAny[massStr].insert(sequence[i->Index()]);
      } else {
        modsVarAny[massStr].insert(sequence[i->Index()]);
      }
      break;
    case PEPTIDE_N:
      if (i->Static()) {
        modsStaticPepN.insert(massStr);
      } else {
        modsVarPepN.insert(massStr);
      }
      break;
    case PEPTIDE_C:
      if (i->Static()) {
        modsStaticPepC.insert(massStr);
      } else {
        modsVarPepC.insert(massStr);
      }
      break;
    case PROTEIN_N:
      carp(CARP_FATAL, "Not yet implemented: PROTEIN_N modifications");
      break;
    case PROTEIN_C:
      carp(CARP_FATAL, "Not yet implemented: PROTEIN_C modifications");
      break;
    }
  }

  double modMass = calcModMass(match);
  string modsSpec = modMass >= 0
    ? "1" + allAminoAcids + "+" + StringUtils::ToString(modMass, modSpecPrecision)
    : "1" + allAminoAcids + StringUtils::ToString(modMass, modSpecPrecision);
  for (map< string, set<char> >::const_iterator i = modsStaticAny.begin(); i != modsStaticAny.end(); i++) {
    modsSpec.push_back(',');
    modsSpec += StringUtils::Join(i->second) + i->first;
  }
  for (map< string, set<char> >::const_iterator i = modsVarAny.begin(); i != modsVarAny.end(); i++) {
    modsSpec.push_back(',');
    modsSpec += StringUtils::ToString(sequence.length()) + StringUtils::Join(i->second) + i->first;
  }

  string ctPep;
  for (set<string>::const_iterator i = modsStaticPepC.begin(); i != modsStaticPepC.end(); i++) {
    if (!ctPep.empty()) {
      ctPep.push_back(',');
    }
    ctPep.push_back(sequence[sequence.length() - 1]);
    ctPep += *i;
  }
  for (set<string>::const_iterator i = modsVarPepC.begin(); i != modsVarPepC.end(); i++) {
    if (!ctPep.empty()) {
      ctPep.push_back(',');
    }
    ctPep.push_back('1');
    ctPep.push_back(sequence[sequence.length() - 1]);
    ctPep += *i;
  }

  string ntPep;
  for (set<string>::const_iterator i = modsStaticPepN.begin(); i != modsStaticPepN.end(); i++) {
    if (!ntPep.empty()) {
      ntPep.push_back(',');
    }
    ntPep.push_back(sequence[0]);
    ntPep += *i;
  }
  for (set<string>::const_iterator i = modsVarPepN.begin(); i != modsVarPepN.end(); i++) {
    if (!ntPep.empty()) {
      ntPep.push_back(',');
    }
    ntPep.push_back('1');
    ntPep.push_back(sequence[0]);
    ntPep += *i;
  }

  VariableModTable* modTable = new VariableModTable();
  if (!modTable->Parse(modsSpec.c_str())) {
    carp(CARP_FATAL, "Error parsing mods");
  } else if (!ctPep.empty() && !modTable->Parse(ctPep.c_str(), CTPEP)) {
    carp(CARP_FATAL, "Error parsing mods (CTPEP)");
  } else if (!ntPep.empty() && !modTable->Parse(ntPep.c_str(), NTPEP)) {
    carp(CARP_FATAL, "Error parsing mods (NTPEP)");
  }
  modTable->SerializeUniqueDeltas();
  return modTable;
}

