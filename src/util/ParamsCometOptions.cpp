#include "Params.h"
#include "AminoAcidUtil.h"
#include "app/CometApplication.h"
#include "app/CometIndexApplication.h"
#include "model/objects.h"
#include "parameter.h"
#include "StringUtils.h"

#include <algorithm>

using namespace std;

void Params::InitCometOptions() {
  /*
   * Comet parameters
   */
  InitArgParam("input spectra",
    "The name of one or more files from which to parse the spectra. Valid formats include mzXML, "
    "mzML, mz5, raw, ms2, and cms2. Files in mzML or mzXML may be compressed with gzip. "
    "RAW files can be parsed only under windows and if the appropriate libraries were "
    "included at compile time. Multiple files can be included on the command line "
    "(space delimited), prior to the name of the database.");
  /* Comet - Database */
  InitArgParam("database_name",
    "A full or relative path to the sequence database to search."
    "The database may be in FASTA or PEFF format, or it may be an index "
    "created with the <code>comet-index</code> command. "
    "Example databases include "
    "RefSeq or UniProt.  The database can contain amino acid "
    "sequences or nucleic acid sequences. If sequences are "
    "amino acid sequences, set the parameter \"nucleotide_reading_frame = 0\". "
    "If the sequences are nucleic acid sequences, you must instruct Comet to "
    "translate these to amino acid sequences. Do this by setting "
    "nucleotide_reading_frame\" to a value between 1 and 9.");
  InitArgParam("index_database_name",
    "A full or relative path to the sequence database to generate "
    "an index for. The database may be in FASTA or PEFF format. "
    "The index will be created in the same directory as the database "
    "and with the same name, but the suffix '.idx'.");
  InitIntParam("decoy_search", 0, 0, 2,
    "0=no, 1=concatenated search, 2=separate search.",
    "Available for comet.", true);
  InitIntParam("peff_format", 0, 0, 5,
    "0=normal FASTA format,"
    "1=PEFF PSI-MOD modifications and amino acid variants,"
    "2=PEFF Unimod modifications and amino acid variants,"
    "3=PEFF PSI-MOD modifications, skipping amino acid variants,"
    "4=PEFF Unimod modifications, skipping amino acid variants,"
    "5=PEFF amino acid variants, skipping PEFF modifications.",
    "Available for comet.", true);
  InitStringParam("peff_obo", "",
    "A full or relative path to the OBO file used with a PEFF search. "
    "Supported OBO formats are PSI-Mod and Unimod OBO files. "
    "Which OBO file you use depends on your PEFF input file. "
    "This parameter is ignored if \"peff_format = 0\". "
    "There is no default value if this parameter is missing.",
    "Available for comet.", true);
  /* Comet - CPU threads */
  InitIntParam("num_threads", 0, -64, 64,
    "0=poll CPU to set num threads; else specify num threads directly.",
    "Available for comet.", true);
  /* Comet - Masses */
  InitDoubleParam("peptide_mass_tolerance", 3.0, 0, BILLION,
    "Controls the mass tolerance value.  The mass tolerance "
    "is set at +/- the specified number i.e. an entered value "
    "of \"1.0\" applies a -1.0 to +1.0 tolerance. "
    "The units of the mass tolerance is controlled by the parameter "
    "\"peptide_mass_units\". ",
    "Available for comet.", true);
  InitDoubleParam("peptide_mass_tolerance_lower", -3.0, -BILLION, BILLION,
    "Controls the lower bound of the precursor mass tolerance value."
    "The units of the mass tolerance is controlled by the parameter "
    "\"peptide_mass_units\". ",
    "Available for comet.", true);
  InitDoubleParam("peptide_mass_tolerance_upper", 3.0, 0, BILLION,
    "Controls the upper bound of the precursor mass tolerance value."
    "The units of the mass tolerance is controlled by the parameter "
    "\"peptide_mass_units\". ",
    "Available for comet.", true);
  InitIntParam("peptide_mass_units", 0, 0, 2,
    "0=amu, 1=mmu, 2=ppm.",
    "Available for comet.", true);
  InitStringParam("auto_peptide_mass_tolerance", "false", "false|warn|fail",
    "Automatically estimate optimal value for the peptide_mass_tolerance parameter "
    "from the spectra themselves. false=no estimation, warn=try to estimate "
    "but use the default value in case of failure, fail=try to estimate and "
    "quit in case of failure.",
    "Available for comet.", true);
  InitIntParam("mass_type_parent", 1, 0, 1,
    "0=average masses, 1=monoisotopic masses.",
    "Available for comet.", true);
  InitIntParam("mass_type_fragment", 1, 0, 1,
    "0=average masses, 1=monoisotopic masses.",
    "Available for comet.", true);
  InitIntParam("precursor_tolerance_type", 0, 0, 1,
    "0=singly charged peptide mass, 1=precursor m/z.",
    "Available for comet.", true);
  InitIntParam("isotope_error", 0, 0, 7,
    "0=off, 1=0/1 (C13 error), 2=0/1/2, 3=0/1/2/3, " 
    "4=-1/0/1/2/3, "
    "5=-1/0/1, "
    "6=-3/-2/-1/0/+1/+2/+3, "
    "7=-8/-4/0/+4/+8 (for +4/+8 stable isotope labeling).",
    "Available for comet.", true);
  /* Comet - Search enzyme */
  InitIntParam("search_enzyme_number", 1, 0, BILLION,
    "Specify a search enzyme from the end of the parameter file.",
    "Available for comet.", true);
  InitIntParam("search_enzyme2_number", 0, 0, BILLION,
    "Specify a second search enzyme from the end of the parameter file.",
    "Available for comet.", true);
  InitIntParam("num_enzyme_termini", 2, 1, 9,
    "valid values are 1 (semi-digested), 2 (fully digested), 8 N-term, 9 C-term.",
    "Available for comet.", true);
  InitIntParam("allowed_missed_cleavage", 2, 0, 5,
    "Maximum value is 5; for enzyme search.",
    "Available for comet.", true);
  /* Comet - Fragment ions */
  InitDoubleParam("fragment_bin_tol", 1.000507, 0, BILLION,
    "Binning to use on fragment ions.",
    "Available for comet.", true);
  InitDoubleParam("fragment_bin_offset", 0.40, 0, 1.0,
    "Offset position to start the binning (0.0 to 1.0).",
    "Available for comet and kojak.", true);
  InitStringParam("auto_fragment_bin_tol", "false", "false|warn|fail",
    "Automatically estimate optimal value for the fragment_bin_tol parameter "
    "from the spectra themselves. false=no estimation, warn=try to estimate "
    "but use the default value in case of failure, fail=try to estimate and "
    "quit in case of failure.",
    "Available for comet.", true);
  InitBoolParam("auto_modifications", false,
    "Automatically infer modifications from the spectra themselves.",
    "Available for comet.", true);
  InitIntParam("theoretical_fragment_ions", 1, 0, 1,
    "0=default peak shape, 1=M peak only.",
    "Available for comet.", true);
  InitIntParam("use_A_ions", 0, 0, 1,
    "Controls whether or not A-ions are considered in the search (0 - no, 1 - yes).",
    "Available for comet.", true);
  InitIntParam("use_B_ions", 1, 0, 1,
    "Controls whether or not B-ions are considered in the search (0 - no, 1 - yes).",
    "Available for comet.", true);
  InitIntParam("use_C_ions", 0, 0, 1,
    "Controls whether or not C-ions are considered in the search (0 - no, 1 - yes).",
    "Available for comet.", true);
  InitIntParam("use_X_ions", 0, 0, 1,
    "Controls whether or not X-ions are considered in the search (0 - no, 1 - yes).",
    "Available for comet.", true);
  InitIntParam("use_Y_ions", 1, 0, 1,
    "Controls whether or not Y-ions are considered in the search (0 - no, 1 - yes).",
    "Available for comet.", true);
  InitIntParam("use_Z_ions", 0, 0, 1,
    "Controls whether or not Z-ions are considered in the search (0 - no, 1 - yes).",
    "Available for comet.", true);
  InitIntParam("use_NL_ions", 1, 0, 1,
    "0=no, 1= yes to consider NH3/H2O neutral loss peak.",
    "Available for comet.", true);
  InitIntParam("use_Z1_ions", 0, 0, 1,
    "Controls whether or not Z1-ions are considered in the search (0 - no, 1 - yes).",
    "Available for comet.", true);
  /* Comet - Output */
  InitIntParam("export_additional_pepxml_scores", 0, 0, 1,
    "0=no, 1=yes Controls whether to output additional search scores in the pep.xml output.",
    "Available for comet.", false);
  InitIntParam("output_mzidentmlfile", 0, 0, 1,
    "0=no, 1=yes  write mzIdentML file.",
    "Available for comet.", true);
  InitIntParam("output_sqtstream", 0, 0, 1,
    "0=no, 1=yes  write sqt file.",
    "Available for comet.", true);
  InitIntParam("output_sqtfile", 0, 0, 1,
    "0=no, 1=yes  write sqt file.",
    "Available for comet.", true);
  InitIntParam("output_txtfile", 1, 0, 1,
    "0=no, 1=yes  write tab-delimited text file.",
    "Available for comet.", true);
  InitIntParam("output_pepxmlfile", 1, 0, 1,
    "0=no, 1=yes  write pep.xml file.",
    "Available for comet.", true);
  InitIntParam("output_percolatorfile", 0, 0, 1,
    "0=no, 1=yes write percolator file.",
     "Available for comet.", true);
  InitIntParam("print_expect_score", 1, 0, 1,
    "0=no, 1=yes to replace Sp with expect in out & sqt.",
    "Available for comet.", false);
  InitIntParam("num_output_lines", 5, 1, BILLION,
    "num peptide results to show.",
    "Available for comet.", true);
  InitIntParam("show_fragment_ions", 0, 0, 1,
    "0=no, 1=yes for out files only.",
    "Available for comet.", false);
  InitIntParam("sample_enzyme_number", 1, 0, 10,
    "Sample enzyme which is possibly different than the one applied to the search. "
    "Used to calculate NTT & NMC in pepXML output.",
    "Available for comet. ", false);
  /* Comet - mzXML/mzML parameters */
  InitStringParam("scan_range", "0 0",
    "Start and scan scan range to search; 0 as first entry ignores parameter.",
    "Available for comet.", true);
  InitStringParam("precursor_charge", "0 0",
    "Precursor charge range to analyze; does not override "
    "mzXML charge; 0 as first entry ignores parameter.",
    "Available for comet.", true);
  InitIntParam("override_charge", 0, 0, 3,
    "Specifies the whether to override existing precursor charge state information when present "
    "in the files with the charge range specified by the \"precursor_charge\" parameter.",
    "Available for comet.", true);
  InitIntParam("ms_level", 2, 2, 3,
    "MS level to analyze, valid are levels 2 or 3.",
    "Available for comet. ", true);
  InitStringParam("activation_method", "ALL", "ALL|CID|ECD|ETD+SA|ETD|PQD|HCD|IRMPD",
    "Specifies which scan types are searched.",
    "Available for comet. ", true);
  /* Comet - Misc. parameters */
  InitStringParam("digest_mass_range", "600.0 5000.0",
    "MH+ peptide mass range to analyze.",
    "Available for comet.", true);
  InitIntParam("num_results", 50, 0, BILLION,
    "Number of search hits to store internally.",
    "Available for comet.", true);
  InitIntParam("skip_researching", 1, 0, 1,
    "For '.out' file output only, 0=search everything again, 1=don't search if .out exists.",
    "Available for comet.", false);
  InitIntParam("max_fragment_charge", 3, 1, 5,
    "Set maximum fragment charge state to analyze (allowed max 5).",
    "Available for comet.", true);
  InitIntParam("max_index_runtime", 0, 0, BILLION,
    "Sets the maximum indexed database search run time for a scan/query. "
    "Valid values are integers 0 or higher representing the maximum run time "
    "in milliseconds. "
    "As Comet loops through analyzing peptides from the database index file, "
    "it checks the cummulative run time of that spectrum search after each "
    "peptide is analyzed. If the run time exceeds the value set for this "
    "parameter, the search is aborted and the best peptide result analyzed "
    "up to that point is returned. "
    "To have no maximum search time, set this parameter value to \"0\". "
    "The default value is \"0\".",
    "Available for comet.", false);
  InitIntParam("max_precursor_charge", 6, 1, 9,
    "Set maximum precursor charge state to analyze (allowed max 9).",
    "Available for comet.", true);
  InitIntParam("nucleotide_reading_frame", 0, 0, 9,
    "0=proteinDB, 1-6, 7=forward three, 8=reverse three, 9=all six.",
    "Available for comet.", true);
  InitIntParam("clip_nterm_methionine", 0, 0, 1,
    "0=leave sequences as-is; 1=also consider sequence w/o N-term methionine.",
    "Available for comet.", true);
  InitIntParam("explicit_deltacn", 0, 0, 1,
    "0=Comet deltaCn reported between the top peptide and the first dissimilar peptide, "
    "1=Comet deltaCn reported between the top two peptides.",
    "Available for comet.", false);
  InitIntParam("old_mods_encoding", 0, 0, 1,
    "0=Comet will use mass based modification encodings, "
    "1=Comet will use the old character based modification encodings.",
    "Available for comet.", false);
  InitIntParam("resolve_fullpaths", 0, 0, 1,
    "Controls whether or not to resolve the full paths of the input files. "
    "0=Comet will not resolve full paths, "
    "1=Comet will resolve full paths.",
    "Available for comet.", true);
  InitStringParam("pinfile_protein_delimiter", "",
    "The default delimiter for the protein field is a tab. "
    "This parameter allows one to specify a different character or string "
    "If this parameter is left blank or is missing, the default tab delimitter is used.",
    "Available for comet.", true);
  InitIntParam("spectrum_batch_size", 20000, 0, BILLION,
    "Maximum number of spectra to search at a time; 0 to search the entire scan range in one loop.",
    "Available for comet.", true);
  InitStringParam("decoy_prefix", "decoy_",
    "Specifies the prefix of the protein names that indicates a decoy.",
    "Available for comet.", true);
  InitStringParam("text_file_extension", "",
    "Specifies the a custom extension for output text file.",
    "Available for comet.", true);
  InitStringParam("output_suffix", "",
    "Specifies the suffix string that is appended to the base output name "
    "for the pep.xml, pin.xml, txt and sqt output files.",
    "Available for comet.", true);
  InitIntParam("peff_verbose_output", 0, 0, 1,
    "Specifies whether the verbose output is reported during a PEFF search. "
    "To show verbose output, set the value to 1. "
    "The default value is 0 if this parameter is missing.",
    "Available for comet.", false);
  InitStringParam("peptide_length_range", "6 50",
    "Defines the length range of peptides to search. "
    "This parameter has two integer values. "
    "The first value is the minimum length cutoff and the second value is "
    "the maximum length cutoff. Only peptides within the specified length "
    "range are analyzed. The maximum peptide length that Comet can analyze is 63. "
    "The default values are \"1 50\".",
    "Available for comet.", true);
  InitStringParam("precursor_NL_ions", "",
    "Controls whether or not precursor neutral loss peaks are considered in "
    "the xcorr scoring. If left blank, this parameter is ignored.  To consider "
    "precursor neutral loss peaks, add one or more neutral loss mass value "
    "separated by a space.  Each entered mass value will be subtracted from "
    "the experimentral precursor mass and resulting neutral loss m/z values "
    "for all charge states (from 1 to precursor charge) will be analyzed. "
    "As these neutral loss peaks are analyzed along side fragment ion peaks, "
    "the fragment tolerance settings (fragment_bin_tol, fragment_bin_offset, "
    "theoretical_fragment_ion) apply to the precursor neutral loss peaks. "
    "The default value is blank/unused.",
    "Available for comet.", true);
  InitIntParam("equal_I_and_L", 1, 0, 1,
    "This parameter controls whether the Comet treats isoleucine (I) and "
    "leucine (L) as the same/equivalent with respect to a peptide identification. "
    "0 treats I and L as different, 1 treats I and L as the same. "
    "The default value is \"1\"",
    "Available for comet.", true);
  InitStringParam("mass_offsets", "",
    "Specifies one or more mass offsets to apply. This value(s) are effectively "
    "subtracted from each precursor mass such that peptides that are smaller "
    "than the precursor mass by the offset value can still be matched to the "
    "respective spectrum.",
    "Available for comet.", true);
  InitIntParam("max_duplicate_proteins", 20, -1, BILLION,
    "defines the maximum number of proteins (identifiers/accessions) to report. "
    "If a peptide is present in 6 total protein sequences, there is one (first) "
    "reference protein and 5 additional duplicate proteins. This parameter "
    "controls how many of those 5 additional duplicate proteins are reported."
    "If \"decoy_search = 2\" is set to report separate target and decoy results, "
    "this parameter will be applied to the target and decoy outputs separately. "
    "If set to \"-1\", there will be no limit on the number of reported additional proteins. "
    "The default value is \"20\" if this parameter is missing.",
    "Available for comet.", true);
  /* Comet - Spectral processing */
  InitIntParam("minimum_peaks", 10, 1, BILLION,
    "Minimum number of peaks in spectrum to search.",
    "Available for comet.", true);
  InitDoubleParam("minimum_intensity", 0, 0, BILLION,
    "Minimum intensity value to read in.",
    "Available for comet. ", true);
  InitIntParam("remove_precursor_peak", 0, 0, 2,
    "0=no, 1=yes, 2=all charge reduced precursor peaks (for ETD).",
    "Available for comet. ", true);
  InitDoubleParam("remove_precursor_tolerance", 1.5, -BILLION, BILLION,
    "+- Da tolerance for precursor removal.",
    "Available for comet. ", true);
  InitStringParam("clear_mz_range", "0.0 0.0",
    "For iTRAQ/TMT type data; will clear out all peaks in the specified m/z range.",
    "Available for comet.", true);
  /* Comet - Variable modifications */
  InitStringParam("variable_mod01", "0.0 null 0 3 -1 0 0",
                  "Up to 15 variable modifications are supported. Each modification "
                  "is specified using seven entries: "
                  "\"[[html:&lt;mass&gt;]][[nohtml:<mass>]] "
                  "[[html:&lt;residues&gt;]][[nohtml:<residues>]] "
                  "[[html:&lt;type&gt;]][[nohtml:<type>]] "
                  "[[html:&lt;max&gt;]][[nohtml:<max>]] "
                  "[[html:&lt;distance&gt;]][[nohtml:<distance>]] "
                  "[[html:&lt;terminus&gt;]][[nohtml:<terminus>]] "
                  "[[html:&lt;force&gt;]][[nohtml:<force>]]\". "
                  "Type is 0 for static mods and non-zero for variable mods. "
                  "Note that that if you set the same type value on multiple "
                  "modification entries, Comet will treat those variable modifications "
                  "as a binary set. This means that all modifiable residues in the "
                  "binary set must be unmodified or modified. Multiple binary sets "
                  "can be specified by setting a different binary modification value. "
                  "Max is an integer specifying the maximum number of modified "
                  "residues possible in a peptide for this modification entry. "
                  "Distance specifies the distance the modification is applied to "
                  "from the respective terminus: -1 = no distance contraint; "
                  "0 = only applies to terminal residue; N = only applies to "
                  "terminal residue through next N residues. "
                  "Terminus specifies which terminus the distance constraint is "
                  "applied to: 0 = protein N-terminus; 1 = protein C-terminus; "
                  "2 = peptide N-terminus; 3 = peptide C-terminus."
                  "Force specifies whether peptides must contain this modification: "
                  "0 = not forced to be present; 1 = modification is required.",
                  "Available for comet.", true);
  for (int i = 2; i <= 9; i++) {
    InitStringParam("variable_mod0" + StringUtils::ToString(i), "0.0 null 0 3 -1 0 0",
                    "See syntax for variable_mod01.",
                    "Available for comet.", true);
  }
  for (int i = 10; i <= 15; i++) {
    InitStringParam("variable_mod" + StringUtils::ToString(i), "0.0 null 0 3 -1 0 0",
                    "See syntax for variable_mod01.",
                    "Available for comet.", true);
  }
  InitIntParam("max_variable_mods_in_peptide", 5, 0, BILLION,
    "Specifies the total/maximum number of residues that can be modified in a peptide.",
    "Available for comet.", true);
  InitIntParam("require_variable_mod", 0, 0, 1,
    "Controls whether the analyzed peptides must contain at least one variable modification.",
    "Available for comet.", true);
  InitStringParam("protein_modlist_file", "",
    "Specify a full or relative path to a protein modifications file. "
    "If this entry points to a modifications file, Comet will parse the modification numbers and protein strings " 
    "from the file and limit the application of the specified variable modifications to the sequence entries that "
    "match the protein string.",
    "Available for comet.", true);
  /* Comet - Static modifications */
  InitDoubleParam("add_Cterm_peptide", 0, 0, BILLION,
    "Specifiy a static modification to the c-terminus of all peptides.",
    "Available for comet.", true);
  InitDoubleParam("add_Nterm_peptide", 0, 0, BILLION,
    "Specify a static modification to the n-terminus of all peptides.",
    "Available for comet.", true);
  InitDoubleParam("add_Cterm_protein", 0, 0, BILLION,
    "Specify a static modification to the c-terminal peptide of each protein.",
    "Available for comet.", true);
  InitDoubleParam("add_Nterm_protein", 0, 0, BILLION,
    "Specify a static modification to the n-terminal peptide of each protein.",
    "Available for comet.", true);
  for (char c = 'A'; c <= 'Z'; c++) {
    InitDoubleParam(CometApplication::staticModParam(c),
                    c != 'C' ? 0 : CYSTEINE_DEFAULT, 
                    -std::numeric_limits<double>::max(), 
                    std::numeric_limits<double>::max(),
                    "Specify a static modification to the residue " + string(1, c) + ".",
                    "Available for comet.", true);
  }
  InitBoolParam("create_peptide_index", false,
    "Create an index of peptides.", "Available for comet-index.", true);
  InitBoolParam("create_fragment_index", false,
    "Create an index of ion fragments.","Available for comet-index",  false);
  InitDoubleParam("fragindex_max_fragmentmass", 2000.0, 0, BILLION,
    "This parameter defines the maximum fragment ion mass to include in the fragment ion index.",
    "Available for comet-index.", true);
  InitDoubleParam("fragindex_min_fragmentmass", 200.0, 0, BILLION,
    "This parameter defines the minimum fragment ion mass to include in the fragment ion index.",
    "Available for comet-index.", true);
  InitIntParam("fragindex_min_ions_report", 3, 1, BILLION,
    "This parameter sets the minimum number fragment ions a peptide must match against the fragment"
    " on index in order to report this peptide in the output",
    "Available for comet-index.", true);
  InitIntParam("fragindex_min_ions_score", 3, 1, BILLION,
    "This parameter sets the minimum number fragment ions a peptide must match against the fragment"
    "ion index in order to proceed to xcorr scoring.",
    "Available for comet-index.", true);
  InitIntParam("fragindex_num_spectrumpeaks", 100, 1, BILLION,
    "This parameter defines the number of mass/intensity pairs that would be queried "
    "against the fragment ion index",
    "Available for comet-index.", true);
  InitIntParam("fragindex_skipreadprecursors", 0, 0, 1,
    "This parameter controls whether or not Comet reads all precursors from the input files. "
    "It uses this information to limit the peptides that are included in the fragment ion index.",
    "Available for comet-index.", true);

  InitBoolParam("list-of-files", false,
    "Specify that the search results are provided as lists of files, rather than as "
    "individual files.",
    "Available for assign-confidence.", true);
}
