#include "LocalizeModification.h"
#include "io/MatchCollectionParser.h"
#include "tide/peptide.h"
#include "tide/modifications.h"
#include "util/FileUtils.h"
#include "util/Params.h"
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
 *     - Convert it to spectrumrecords format
 *     - Load spectra into memory
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

  map<string, string> spectrumFiles; // match file path -> spectrumrecords file
  map<string, string> spectrumFilesRev; // spectrumrecords file -> match file path
  set<string> tempSr; // temp converted spectrumrecords files
  map<string, SpectrumCollection*> preloadedSpectra; // match file path -> SpectrumCollection*

  int numPsms = 0;
  bool hasTargets = false;
  bool hasDecoys = false;

  MatchIterator* matchIter = new MatchIterator(matches);
  while (matchIter->hasNext()) {
    ++numPsms;
    Crux::Match* match = matchIter->next();
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
    if (!SpectrumRecordSpectrumCollection::IsSpectrumRecordFile(spectrumFile)) {
      // convert to spectrumrecords
      carp(CARP_INFO, "Converting %s to spectrumrecords format", matchPath.c_str());
      string convertedFile = FileUtils::Join(
        Params::GetString("output-dir"),
        getName() + ".spectrumrecords" + StringUtils::ToString(tempSr.size()) + ".tmp");
      if (!SpectrumRecordWriter::convert(spectrumFile, convertedFile)) {
        carp(CARP_FATAL, "Error converting %s to spectrumrecords format", spectrumFile.c_str());
      }
      spectrumFile = convertedFile;
      tempSr.insert(convertedFile);
    }
    spectrumFiles[matchPath] = spectrumFile;
    spectrumFilesRev[spectrumFile] = matchPath;
    preloadedSpectra[spectrumFile] = TideSearchApplication::loadSpectra(spectrumFile);
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
  ofstream* outfile = create_stream_in_path(outpath.c_str(), NULL, Params::GetBool("overwrite"));

  carp(CARP_INFO, "Scoring modified peptides (results will be written to %s)...",
       FileUtils::BaseName(outpath).c_str());
  TideMatchSet::writeHeaders(outfile, false, Params::GetBool("compute-sp"));

  int curPsm = 0;
  matchIter = new MatchIterator(matches);
  while (matchIter->hasNext()) {
    Crux::Match* match = matchIter->next();

    // create temp index
    string tempIndex = FileUtils::Join(Params::GetString("output-dir"), getName() + ".tempindex");
    FileUtils::Mkdir(tempIndex);
    if (!FileUtils::IsDir(tempIndex)) {
      carp(CARP_FATAL, "Couldn't create temp directory %s", tempIndex.c_str());
    }
    string protix = FileUtils::Join(tempIndex, "protix");
    string pepix = FileUtils::Join(tempIndex, "pepix");
    string auxlocs = FileUtils::Join(tempIndex, "auxlocs");
    writeIndexProteins(protix, match->getPeptide());
    writeIndexPeptides(pepix, auxlocs, match);

    // run tide-search
    int scan = match->getSpectrum()->getFirstScan();
    carp(CARP_DETAILED_INFO, "Scoring modified forms of %s against spectrum %d",
         match->getPeptide()->getModifiedSequenceWithMasses().c_str(), scan);
    TideSearchApplication search;
    search.localizeMod_ = true;
    search.targetFile_ = outfile;
    search.scanNumber_ = StringUtils::ToString(scan);
    search.spectrumFilesOverride_ = &spectrumFilesRev;
    search.spectra_ = preloadedSpectra;
    search.main(vector<string>(1, spectrumFiles[match->getFilePath()]), tempIndex);

    // remove temp index
    FileUtils::Remove(tempIndex);

    reportProgress(++curPsm, numPsms);
  }
  delete matchIter;
  delete matches;

  delete outfile;

  for (set<string>::const_iterator i = tempSr.begin(); i != tempSr.end(); i++) {
    FileUtils::Remove(*i);
  }

  for (map<string, SpectrumCollection*>::const_iterator i = preloadedSpectra.begin();
       i != preloadedSpectra.end();
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

bool LocalizeModificationApplication::hidden() const { return false; }

void LocalizeModificationApplication::reportProgress(int curPsm, int numPsms) {
  int percent = (int)((double)curPsm / (double)numPsms * 100.0);
  set<int>::iterator i = progress_.find(percent);
  if (i != progress_.end()) {
    progress_.erase(i);
    carp(CARP_INFO, "%d%% complete", percent);
  }
}

void LocalizeModificationApplication::writeIndexProteins(
  const string& proteinsFile,
  Crux::Peptide* peptide
) const {
  pb::Header header;
  header.set_file_type(pb::Header::RAW_PROTEINS);
  HeadedRecordWriter writer(proteinsFile, header);
  int proteinId = -1;
  for (PeptideSrcIterator i = peptide->getPeptideSrcBegin();
       i != peptide->getPeptideSrcEnd();
       i++) {
    pb::Protein protein;
    protein.set_id(++proteinId);
    Crux::Protein* cruxProtein = (*i)->getParentProtein();
    protein.set_name(cruxProtein->getId());
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
    protein.set_residues(residues);
    writer.Write(&protein);
  }
}

void LocalizeModificationApplication::writeIndexPeptides(
  const string& peptidesFile,
  const string& auxLocsFile,
  Crux::Match* match
) const {
  Crux::Peptide* cruxPeptide = match->getPeptide();
  VariableModTable* modTable = getModTable(match);

  pb::Header header;
  header.set_file_type(pb::Header::PEPTIDES);
  pb::Header_PeptidesHeader& pepHeader = *(header.mutable_peptides_header());
  pepHeader.set_min_mass(match->getNeutralMass());
  pepHeader.set_max_mass(match->getNeutralMass());
  pepHeader.set_min_length(cruxPeptide->getLength());
  pepHeader.set_max_length(cruxPeptide->getLength());
  //pepHeader.set_enzyme();
  //pepHeader.set_full_digestion();
  // TODO Detect mono/average from results file?
  pepHeader.set_monoisotopic_precursor(Params::GetString("isotopic-mass") == "mono");
  pepHeader.set_has_peaks(false);
  pepHeader.mutable_mods()->CopyFrom(*(modTable->ParsedModTable()));
  pepHeader.mutable_nterm_mods()->CopyFrom(*(modTable->ParsedNtpepModTable()));
  pepHeader.mutable_cterm_mods()->CopyFrom(*(modTable->ParsedCtpepModTable()));
  //pepHeader.mutable_peptides_header()->set_decoys()
  HeadedRecordWriter writer(peptidesFile, header);

  pb::Header auxLocsHeader;
  auxLocsHeader.set_file_type(pb::Header::AUX_LOCATIONS);
  pb::Header_Source* source = auxLocsHeader.add_source();
  source->set_filename(peptidesFile);
  source->mutable_header()->CopyFrom(header);

  HeadedRecordWriter writer2(auxLocsFile, auxLocsHeader);
  int auxLocationsIndex = -1;

  pb::Peptide peptide;
  peptide.set_mass(match->getNeutralMass());
  peptide.set_length(cruxPeptide->getLength());
  peptide.set_is_decoy(false);

  char* seq = cruxPeptide->getSequence();
  string sequence(seq);
  free(seq);
  
  // add existing modifications
  vector<Crux::Modification> existingMods = cruxPeptide->getVarMods();
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
      if (abs(i->DeltaMass() - modTable->PossDelta(modIdx)) < 0.1) {
        varModIdx = modIdx;
        break;
      }
    }
    if (varModIdx == -1) {
      carp(CARP_FATAL, "Couldn't find mod %.1f in VariableModTable", i->DeltaMass());
    }
    peptide.add_modifications(modTable->EncodeMod((int)i->Index(), varModIdx));
  }

  double modMass = calcModMass(match);
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

  set<unsigned char> existingIndices;
  for (vector<Crux::Modification>::const_iterator i = existingMods.begin();
       i != existingMods.end();
       i++) {
    existingIndices.insert(i->Index());
  }

  int peptideId = -1;
  for (unsigned char i = 0; i < cruxPeptide->getLength(); i++) {
    if (existingIndices.find(i) != existingIndices.end()) {
      continue; // don't modify same AA twice
    }
    pb::AuxLocation locations;
    peptide.set_id(++peptideId);
    peptide.add_modifications(modTable->EncodeMod((int)i, openModIdx));

    int proteinId = -1;
    for (PeptideSrcIterator i = cruxPeptide->getPeptideSrcBegin();
         i != cruxPeptide->getPeptideSrcEnd();
         i++) {
      int start = (*i)->getStartIdxOriginal();
      if (++proteinId == 0) {
        // first src
        peptide.mutable_first_location()->set_protein_id(proteinId);
        peptide.mutable_first_location()->set_pos(start);
      } else {
        // aux location
        pb::Location* location = locations.add_location();
        location->set_protein_id(proteinId);
        location->set_pos(start);
      }
    }
    if (locations.location_size() > 0) {
      peptide.set_aux_locations_index(++auxLocationsIndex);
      writer2.Write(&locations);
    }
    writer.Write(&peptide);
    peptide.mutable_modifications()->RemoveLast();
  }
  delete modTable;
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

