#include "AminoAcidUtil.h"
#include "mass.h"
#include "model/Peptide.h"
#include "objects.h"
#include "parameter.h"
#include "Params.h"

#include <algorithm>

using namespace std;

Params::ParamContainer Params::container_;

void Params::Initialize() {
  if (!container_.Empty() || container_.Finalized()) {
    throw runtime_error("Parameters already initialized");
  }
  /* generate_peptide arguments */
  InitArgParam("protein fasta file",
   "File containing protein sequences in fasta format.");
  InitArgParam("index name",
    "Name to give the new directory containing index files.");
  InitArgParam("ms2 file",
    "File containing spectra to be searched.");
  /* get-ms2-spectrum */
  InitIntParam("scan number", 0, 0, BILLION,
    "Scan number identifying the spectrum.",
    "Argument for get-ms2-spectrum", false);
  InitArgParam("output file",
    "File where spectrum will be written.");
  /* predict-peptide-ions */
  InitArgParam("peptide sequence",
    "The sequence of the peptide.");
  InitArgParam("charge state",
    "The charge state of the peptide.");
  /* hardklor arguments */
  InitArgParam("spectra",
    "The name of a file from which to parse high-resolution spectra. The file "
    "may be in MS1 (.ms1), binary MS1 (.bms1), compressed MS1 (.cms1), or mzXML "
    "(.mzXML) format.");
  /*Percolator arguments*/
  InitArgParam("pin",
    "PIN files are tab-delimited files for PIN format. "
    "Also, this argument can be \"-\" which indicates the pin file will come "
    "from standard input. Alternately, a SQT, PepXML, or tab-delimited file may "
    "be given (a corresponding decoy file must also exist in the same directory), "
    "in which case a pin file will be generated in the output directory prior "
    "to execution.");
  /*make-pin arguments*/
  InitArgParam("target input",
    "search results file in sqt, tab-delimited or pep.xml format.  "
    "Also, this argument can be - which indicates the result file will come from "
    "standard input");
  InitStringParam("decoy input", "",
    "make-pin can convert any file format in sqt, tab-delimited and pep.xml file "
    "to pin file ",
    "Argument, not option for make-pin", false);
  InitStringParam("output-file", "",
    "Path where pin file will be written",
    "It is optional for make-pin", false);
  InitBoolParam("filestem-prefixes", false,
    "Prefix PSM IDs with filestems instead of target or decoy and file index.",
    "Available for make-pin", false);
  InitBoolParam("mod-symbols", false,
    "Print modification symbols instead of masses in peptide sequences.",
    "Available for make-pin", false);
  /* *** Initialize Options (command line and param file) *** */

  /* options for all executables */
  InitBoolParam("version", false, "Print version number and quit.",
    "Available for all crux programs.  On command line use '--version T'.", true);
  InitIntParam("verbosity", 30, 0, 100,
    "Set level of output to stderr (0-100).  Default=30.",
    "Available for all crux programs.  Each level prints the following "
    "messages, including all those at lower verbosity levels: 0-fatal "
    "errors, 10-non-fatal errors, 20-warnings, 30-information on the "
    "progress of execution, 40-more progress information, 50-debug info, "
    "60-detailed debug info.", true);
  InitStringParam("parameter-file", "", 
    "Set additional options with values in the given file.",
    "Available for all crux programs. Any options specified on the "
    "command line will override values in the parameter file.", true);
  InitBoolParam("overwrite", false, 
    "Replace existing files (T) or exit if attempting to "
    "overwrite (F). Default=F.",
    "Available for all crux programs.  Applies to parameter file "
    "as well as index, search, and analysis output files.", true);
  /* generate_peptide parameters  */
  InitIntParam("min-length", 6, 1, MAX_PEPTIDE_LENGTH,
    "The minimum length of peptides to consider. Default=6.",
    "Used from the command line or parameter file by "
    "crux-generate-peptides, crux tide-index, and crux generate-decoys.", true);
  InitIntParam("max-length", 50, 1, MAX_PEPTIDE_LENGTH,
    "The maximum length of peptides to consider. Default=50.",
    "Available from command line or parameter file for "
    "crux-generate-peptides, crux tide-index, and crux generate-decoys. ", true);
  InitDoubleParam("min-mass", 200, 0, BILLION,
    "The minimum mass of peptides to consider. Default=200.",
    "Available from command line or parameter file for "
    "crux-generate-peptides, crux tide-index, and crux generate-decoys. ", true);
  InitDoubleParam("max-mass", 7200, 1, BILLION, 
    "The maximum mass of peptides to consider. Default=7200.",
    "Available from command line or parameter file for "
    "crux-generate-peptides, crux tide-index, and crux generate-decoys. ", true);
  InitStringParam("isotopic-mass", "average", "average|mono",
    "Which isotopes to use in calcuating peptide mass. "
    "<string>=average|mono. Default=average.", 
    "Used from command line or parameter file by "
    "crux-generate-peptides and crux generate-decoys.", true);
  InitIntParam("min-peaks", 20, 0, BILLION,
    "The minimum number of peaks a spectrum must have for it to be searched."
    " Default=20.", 
    "Available for tide-search.", true);
  InitStringParam("digestion", "full-digest",
    "full-digest|partial-digest|non-specific-digest",
    "Degree of digestion used to generate peptides. "
    "<string>=full-digest|partial-digest. Either both ends or one end "
    "of a peptide must conform to enzyme specificity rules. "
    "Default=full-digest.",
    "Used in conjunction with enzyme option when enzyme is not set to "
    "to 'no-enzyme'.  Available from command line or parameter file for "
    "crux-generate-peptides, crux tide-index, and crux "
    "generate-decoys. Digestion rules are as "
    "follows: enzyme name [cuts after one of these residues][but not before "
    "one of these residues].  trypsin [RK][P], elastase [ALIV][P], "
    "chymotrypsin [FWYL][P].", true);
  InitStringParam("enzyme", "trypsin", "no-enzyme|trypsin|trypsin/p|chymotrypsin|"
    "elastase|clostripain|cyanogen-bromide|iodosobenzoate|proline-endopeptidase|"
    "staph-protease|asp-n|lys-c|lys-n|arg-c|glu-c|pepsin-a|"
    "elastase-trypsin-chymotrypsin|custom-enzyme",
    "Enzyme to use for in silico digestion of proteins. "
    "<string>=no-enzyme|trypsin|trypsin/p|chymotrypsin| " 
    "elastase|clostripain|cyanogen-bromide|iodosobenzoate| " 
    "proline-endopeptidase|staph-protease|asp-n|lys-c "
    "lys-n|arg-c|glu-c|pepsin-a| "
    "|elastase-trypsin-chymotrypsin|custom-enzyme. "
    "Default=trypsin.", 
    "Used in conjunction with the options digestion and missed-cleavages. "
    "Use 'no-enzyme' for non-specific digestion.  Available "
    "from command line or parameter file for crux-generate-peptides, "
    "crux tide-index, and crux generate-decoys.  "
    "Digestion rules: enzyme name [cuts after one of these residues]|{but "
    "not before one of these residues}. trypsin [RK]|{P}, trypsin/p [RK]|[], "
    "elastase [ALIV]|{P}, chymotrypsin [FWYL]|{P}, clostripain [R]|[], "
    "cyanogen-bromide [M]|[], iodosobenzoate [W]|[], proline-endopeptidase "
    "[P]|[], staph-protease [E]|[], elastase-trypsin-chymotrypsin "
    "[ALIVKRWFY]|{P},asp-n []|[D] (cuts before D), lys-c [K]|{P}, lys-n "
    "[]|[K] (cuts before K), arg-c [R]|{P}, glu-c [DE]|{P}, pepsin-a "
    "[FL]|{P}.", true);
  InitStringParam("precursor-window-type", "mass", "mass|mz|ppm",
    "Window type to use for selecting candidate "
    "peptides.  <string>=mass|mz|ppm. Default=mass.",
    "Available for search-for-xlinks and tide-search.", true);
  InitStringParam("spectrum-parser", "pwiz", "pwiz|mstoolkit",
    "Parser to use for reading in spectra "
    "<string>=pwiz|mstoolkit. Default=pwiz.",
    "Available for search-for-xlinks.", true);
  InitStringParam("custom-enzyme", "", 
    "Specify rules for in silico digestion of proteins. "
    "See HTML documentation for syntax. Default is trypsin.",
    "Overrides the enzyme option.  Two lists of residues are given enclosed "
    "in square brackets or curly braces and separated by a |. The first list "
    "contains residues required/prohibited before the cleavage site and the "
    "second list is residues after the cleavage site.  If the residues are "
    "required for digestion, they are in square brackets, '[' and ']'.  "
    "If the residues prevent digestion, then they are enclosed in curly "
    "braces, '{' and '}'.  Use X to indicate all residues.  For example, "
    "trypsin cuts after R or K but not before P which is represented as "
    "[RK]|{P}.  AspN cuts after any residue but only before D which is "
    "represented as [X]|[D].", true);
  InitIntParam("missed-cleavages", 0, 0, 500,
    "Include peptides with up to n missed cleavage sites. Default=0.",
    "Available from command line or parameter file for "
    "crux-generate-peptides and crux generate-decoys. "
    "When used with enzyme=<trypsin|elastase|chymotrypsin> "
    "includes peptides containing one or more potential cleavage sites.", true);
  InitStringParam("keep-terminal-aminos", "NC", 
    "When creating decoy peptides using decoy-format=shuffle or decoy-format="
    "peptide-reverse, this option specifies whether the N-terminal and "
    "C-terminal amino acids are kept in place or allowed to be shuffled or "
    "reversed. Default = NC.",
    "Available for tide-index and generate-decoys.", true);
  InitBoolParam("unique-peptides", true,
    "Generate peptides only once, even if they appear in more "
    "than one protein (T,F).  Default=T.",
    "Available from command line or parameter file for "
    "crux-genereate-peptides. Returns one line per peptide "
    "when true or one line per peptide per protein occurence when false.  ",
    true);
  InitBoolParam("peptide-list", false,
    "Create an ASCII version of the peptide list. Default=F.",
    "Creates an ASCII file in the output directory "
    "containing one peptide per line.", true);
   // print-processed-spectra option
  InitStringParam("stop-after", "xcorr", "remove-precursor|square-root|"
    "remove-grass|ten-bin|xcorr",
    "Specify which preprocessing step to stop after. Value must be discretize, "
    "remove-precursor, square-root, remove-grass, ten-bin, or xcorr. Default = "
    "xcorr.",
    "Available for print-processed-spectra.", true);
  /* more generate_peptide parameters */
  InitBoolParam("output-sequence", false, 
    "Print peptide sequence (T,F). Default=F.",
    "Available only for crux-generate-peptides.", true);
  InitBoolParam("sqt-output", false,
    "Output SQT in the output directory.  Default=F",
    "Available for tide-search.", true);
  InitBoolParam("mzid-output", false,
    "Output MZID in the output directory.  Default=F",
    "Available for tide-search, percolator.", true);
  InitBoolParam("pin-output", false,
    "Output PIN XML in the output directory.  Default=F",
    "Available for tide-search.", true);
  InitBoolParam("pout-output", false,
    "Output POUT XML in the output directory.  Default=F",
    "Available for percolator.", true);
  InitBoolParam("pepxml-output", false,
    "Output pepXML in the output directory.  Default=F",
    "Available for tide-search, q-ranker, barista, percolator.", true);
  InitBoolParam("txt-output", true,
    "Output tab-delimited text in the output directory.  Default=T",
    "Available for tide-saerch, percolator, q-ranker, barista.", true);
  InitStringParam("prelim-score-type", "sp", "sp|xcorr",
    "Initial scoring (sp, xcorr). Default=sp,", 
    "The score applied to all possible psms for a given spectrum. Typically "
    "used to filter out the most plausible for further scoring. See "
    "max-rank-preliminary and score-type.", false);
  InitStringParam("score-type", "xcorr", "xcorr|sp|xcorr-pvalue|sp-pvalue",
    "The primary scoring method to use (xcorr, sp, xcorr-pvalue, sp-pvalue)."
    " Default=xcorr.", 
    "Primary scoring is "
    "typically done on a subset (see max-rank-preliminary) of all "
    "possible psms for each spectrum. Default is the SEQUEST-style xcorr."
    " Crux also offers a p-value calculation for each psm based on xcorr "
    "or sp (xcorr-pvalue, sp-pvalue).", false);
  InitBoolParam("compute-sp", false,
    "Compute the Sp score for all candidate peptides.  Default=F",
    "Available for tide-search.", true);
  InitBoolParam("compute-p-values", false, 
    "Compute p-values for the main score type. Default=F.",
    "Currently only implemented for XCORR.", true);
  InitStringParam("scan-number", "",
    "Search only select spectra specified as a single "
    "scan number or as a range as in x-y.  Default=search all.",
    "The search range x-y is inclusive of x and y.", true);
  /* N.B. Use NaN to indicate that no user preference was specified.
   * In this case, the default value depends on the mass type.
   * S.M. Also prevent a width of 0.                                */
  InitDoubleParam("mz-bin-width", 1.0005079, 1e-4, BILLION,
    "Specify the width of the bins used to "
    "discretize the m/z axis.  Also used as tolerance for assigning "
    "ions.  Default=1.0005079 for monoisotopic mass "
    "or 1.0011413 for average mass.",
    "Available for tide-search and xlink-assign-ions.", true);
  InitDoubleParam("mz-bin-offset", 0.40, 0.0, 1.0,
    "Specify the location of the left edge of the "
    "first bin used to discretize the m/z axis. Default=0.40",
    "Available for tide-search.", true);
  // initialize as "unset", then set as bool after cmdline parsed
  InitBoolParam("use-flanking-peaks", false,
    "Include peaks +/- 1da around b/y ions in theoretical spectrum. default=F.",
    "Available for the tide-search and search-for-xlinks commands.", true);
  InitDoubleParam("spectrum-min-mz", 0.0, 0, BILLION, 
    "The lowest spectrum m/z to search. Default=0.0.",
    "Available for tide-search.", true);
  InitDoubleParam("spectrum-max-mz", BILLION, 1, BILLION, 
    "The highest spectrum m/z to search. Default=no maximum.",
    "Available for tide-search.", true);
  InitStringParam("spectrum-charge", "all", "1|2|3|all",
    "Spectrum charge states to search. <string>=1|2|3|all. Default=all.",
    "Used by tide-search to limit the charge "
    "states considered in the search.  With 'all' every spectrum will be "
    "searched and spectra with multiple charge states will be searched "
    "once at each charge state.  With 1, 2 ,or 3 only spectra with that "
    "that charge will be searched.", true);
  InitStringParam("fileroot", "", 
    "Prefix added to output file names. Default=none. ",
    "Used by crux percolator, crux compute-q-values, and crux q-ranker.", true);
  InitStringParam("output-dir", "crux-output",
    "Folder to which results will be written. Default='crux-output'.",
    "Used by crux compute-q-values, and crux percolator.", true);
  // user options regarding decoys
  InitStringParam("decoys", "peptide-shuffle",
    "Include a decoy version of every peptide by shuffling or reversing the "
    "target sequence.  <string>=none|reverse|protein-shuffle|peptide-shuffle."
    " Use 'none' for no decoys.  Default=peptide-shuffle.",
    "", true);
  InitIntParam("num-decoys-per-target", 1, 0, 10,
    "Number of decoy peptides to search for every target peptide searched."
    "Only valid for fasta searches when --decoys is not none. Default=1.",
    "Use --decoy-location to control where they are returned (which "
    "file(s)) and --decoys to control how targets are randomized.", true);
  InitStringParam("decoy-location", "separate-decoy-files",
    "Specify location of decoy search results. "
    "<string>=target-file|one-decoy-file|separate-decoy-files. "
    "Default=separate-decoy-files.",
    "Applies when decoys is not none.  Use 'target-file' to mix "
    "target and decoy search results in one file. 'one-decoy-file' will "
    "return target results in one file and all decoys in another. "
    "'separate-decoy-files' will create as many decoy files as "
    "num-decoys-per-target.", true);
  // coder options regarding decoys
  InitIntParam("num-decoy-files", 2, 0, 10,
    "Replaces number-decoy-set.  Determined by decoy-location"
    " and num-decoys-per-target",
    "", false);
  InitBoolParam("tdc", false,
    "Target-decoy competition. puts decoy psms in target file. ",
    "Now hidden from the user", false);
  InitBoolParam("decoy-p-values", false,
    "Store all decoy p-values in a file",
    "", false);
  InitIntParam("max-rank-preliminary", 500, 0, BILLION, 
    "Number of psms per spectrum to score with xcorr after preliminary "
    "scoring with Sp. "
    "Set to 0 to score all psms with xcorr. Default=500.",
    "For positive values, the Sp "
    "(preliminary) score acts as a filter; only high scoring psms go "
    "on to be scored with xcorr.  This saves some time.  If set to 0, "
    "all psms are scored with both scores. ", true);
  InitIntParam("top-match", 5, 1, BILLION, 
    "The number of PSMs per spectrum writen to the output " 
    " file(s).  Default=5.",
    "Available for tide-search and crux percolator",
    true);
  InitStringParam("seed", "1",
    "When given a unsigned integer value seeds the random number generator with that value. "
    "When given the string \"time\" seeds the random number generator with the system time. "
    "Default = 1.",
    "Available for all Crux commands.", true);
  InitDoubleParam("precursor-window", 3.0, 0, 100, 
    "Search peptides within +/- 'precursor-window' "
    "of the spectrum mass.  Definition of precursor window depends "
    "upon precursor-window-type. Default=3.0.",
    "Available for tide-search and crux-generate-peptides.", true);
  InitStringParam("fragment-mass", "mono", "average|mono",
    "Which isotopes to use in calculating fragment ion mass. "
    "<string>=average|mono. Default=mono.", 
    "Used by crux-predict-peptide-ions.", true);
  InitStringParam("mod", "NO MODS",
    "Specify a variable modification to apply to peptides.  " 
    "<mass change>:<aa list>:<max per peptide>:<prevents cleavage>:<prevents cross-link>."
    "  Sub-parameters prevents cleavage and prevents cross-link are optional (T|F)."
    "Default=no mods.",
    "Available from parameter file for crux-generate-peptides and "
    "the same must be used for crux compute-q-value.", true);
  InitStringParam("cmod", "NO MODS",
    "Specify a variable modification to apply to C-terminus of peptides. " 
    "<mass change>:<max distance from protein c-term (-1 for no max)>. " 
    "Default=no mods.",       
    "Available from parameter file for crux-generate-peptides and "
    "the same must be used for crux compute-q-value.", true);
  InitStringParam("nmod", "NO MODS",
    "Specify a variable modification to apply to N-terminus of peptides.  " 
    "<mass change>:<max distance from protein n-term (-1 for no max)>",
    "Available from parameter file for crux-generate-peptides and "
    "the same must be used for crux compute-q-value.", true);
  InitIntParam("min-mods", 0, 0, MAX_PEPTIDE_LENGTH,
    "The minimum number of modifications that can be applied to a single " 
    "peptide.  Default=0.",
    "Available for tide-index.", true);
  InitIntParam("max-mods", MAX_PEPTIDE_LENGTH, 0, MAX_PEPTIDE_LENGTH,
    "The maximum number of modifications that can be applied to a single " 
    "peptide.  Default=no limit.",
    "Available for tide-index.", true);
  InitIntParam("max-aas-modified", MAX_PEPTIDE_LENGTH, 0, MAX_PEPTIDE_LENGTH,
    "The maximum number of modified amino acids that can appear in one "
    "peptide.  Each aa can be modified multiple times.  Default=no limit.",
    "", true);
  InitStringParam("mod-mass-format", "mod-only", "mod-only|total|separate",
    "Print the masses of modified sequences in one of three ways 'mod-only', "
    "'total' (residue mass plus modification), or 'separate' (for multiple "
    "mods to one residue): Default 'mod-only'.",
    "Available for generate-peptides.", true);
  InitIntParam("mod-precision", 2, 0, 20,//arbitrary
    "Set the precision for modifications as written to .txt files.",
    "Also changes mods written to parameter file. Set internally based on "
    "the max mod precision in the param file.",
    false);
  InitIntParam("precision", 8, 1, 100, //max is arbitrary
    "Set the precision for scores written to sqt and text files. Default=8.",
    "Available from parameter file for percolator and compute-q-values.", true);
  InitIntParam("mass-precision", 4, 1, 100, // max is arbitrary
    "Set the precision for masses and m/z written to sqt and .txt files. Default=4",
    "Available from parameter file for all commands.", true);
  InitIntParam("print-search-progress", 1000, 0, BILLION,
    "Show search progress by printing every n spectra searched.  Default="
    "1000.", "Set to 0 to show no search progress.  Available for tide-search",
    true);
  // Sp scoring params
  InitDoubleParam("max-mz", 4000, 0, BILLION, 
    "Used in scoring sp.",
    "Hide from users", false);
  InitDoubleParam("fraction-top-scores-to-fit", 0.55, 0, 1, 
    "The fraction of psms per spectrum to use for estimating the "
    "score distribution for calculating p-values. "
    "Not compatible with 'number-top-scores-to-fig'. Default=0.55.",
    "For developers/research only.", false);
  /* analyze-matches options */
  InitStringParam("algorithm", "percolator", "percolator|curve-fit|none",
    "The analysis algorithm to use (percolator, curve-fit, none)."
    " Default=percolator.",
    "Available only for crux-analyze-matches.  Using 'percolator' will "
    "assign a q-value to the top-ranking psm for each spectrum based on "
    "the decoy searches.  Using 'curve-fit' will assign a q-value to same "
    "using the p-values calculated with score-type=<xcorr-pvalue|"
    "sq-pvalue>.  Incorrect combinations of score-type and algorithm cause"
    " undefined behavior. Using 'none' will turn the binary .csm files "
    "into text.", false);
  // **** percolator options. ****
  InitBoolParam("feature-file", false,
    "Optional file into which psm features are printed. Default=F.",
    "Available for percolator and q-ranker.  File will be named "
    "<fileroot>.percolator.features.txt or <fileroot>.qranker.features.txt.", true);
  InitBoolParam("protein", false,
    "output protein level probability. Default=F",
    "Available for crux percolator", true);
  InitBoolParam("decoy-xml-output", false,
    "Include decoys (PSMs, peptides, and/or proteins) in the "
    "xml-output. Only available if -X is used. Default=F",
    "Available for crux percolator", true);
  InitStringParam("decoy-prefix", "decoy_",
    "Option for single SQT file mode defining the name pattern "
    "used for decoy database. Default=decoy_.",
    "Available for percolator", true);
  InitDoubleParam("c-pos", 0.01,
    "Penalty for mistakes made on positive examples. Set by "
    "cross validation if not specified. Default=cross-validate ",
    "Available for crux percolator", true);
  InitDoubleParam("c-neg", 0.0, 0.0, 0.90,
    "Penalty for mistake made on negative examples. Set by cross "
    "validation if not specified or --c-pos not specified.",
    "Available for crux percolator", true);
  InitDoubleParam("test-fdr", 0.01, 0.0, 1.0,
    "False discovery rate threshold for evaluating best cross validation result "
    "and the reported end result. Default is 0.01.",
    "Availble for crux percolator.", true);
  InitIntParam("maxiter", 10, 0, 100000000,
    "Maximum number of iterations for training (default 10).",
    "Available for crux percolator", false);
  InitDoubleParam("train-ratio", 0.6, 0.0, 1.0,
    "Fraction of the negative data set to be used as train set when only providing"
    " one negative set, remaining examples will be used as test set.Default 0.6",
    "Available for crux percolator.", true);
  InitBoolParam("output-weights", false,
    "Output final weights to percolator.target.weights.txt.Default=T.",
    "Available for crux percolator", true);
  InitStringParam("input-weights", "",
    "Read initial weights from the given file (one per line). Default=F",
    "Available for crux percolator ", true);
  InitStringParam("default-direction", "",
    "The most informative feature given as the feature name.  The name can be "
    "preceded by a hyphen (e.g., \"-XCorr\") to indicate that a lower value is better."
    " By default, Percolator will select the feature that produces"
    " the largest set of target PSMs at a specified FDR threshold"
    " (cf. --train-fdr).",
    "Available for crux percolator", true);
  InitBoolParam("unitnorm", false,
    "Use unit normalization [0-1] instead of standard deviation normalization",
    "Available for crux percolator.", true);
  InitDoubleParam("train-fdr", 0.01, 0, BILLION,
    "False discovery rate threshold to define positive examples in training. "
    "Set by cross validation if 0. Default is 0.01.",
    "Available for crux percolator",
    true);
  InitDoubleParam("alpha", 0.0, 0.0, 1.0,
    "Probability with which a present protein emits an associated peptide (--protein T "
    "must be set). Set by grid search if not specified.",
    "Available for crux percolator if --protein T is set.", true);
  InitDoubleParam("beta", 0.0, 0.0, 10.0,
    "Probability of the creation of a peptide from noise (--protein T "
    "must be set). Set by grid search if not specified.",
    "Available for crux percolator if --protein T is set.", true);
  InitDoubleParam("gamma", 0.0, 0.0, 10.0,
    "Prior probability of that a protein is present in the sample (--protein T "
    "must be set). Set by grid search if not specified.",
    "Available for crux percolator if --protein T is set.", true);
  InitBoolParam("test-each-iteration", false,
    "Measure performance on test set each iteration",
    "Available for crux percolator.", true);
  InitBoolParam("static-override", false,
    "Override error check and do not fall back on default score vector in case "
    "of suspect score vector.",
    "Available for crux percolator.", true);
  InitBoolParam("klammer", false,
    "Using retention time features calculated as in Klammer et al.",
    "Available for crux percolator", true);
  InitIntParam("doc", -1, -1, 15,
    "Include description of correct features.",
    "Avilable for crux percolator", true);
  InitBoolParam("only-psms", false,
    "Do not remove redundant peptides, keep all PSMs and exclude peptide level probability.",
    "Available for crux percolator", true);
  InitBoolParam("allow-protein-group", false,
    "Treat ties as if it were one protein ",
    "Available for crux percolator.", true);
  InitBoolParam("protein-level-pi0", false,
    "Use pi_0 value when calculating empirical q-values (--protein T must be set).",
    "Available for crux percolator if --protein T is set.", true);
  InitBoolParam("group-proteins", false,
    "Proteins with same probabilities will be grouped (--protein T must be set).",
    "Available for crux percolator if --protein T is set.", true);
  InitBoolParam("empirical-protein-q", false,
    "Output empirical q-values from target-decoy analysis (--protein T must be set).",
    "Available for crux percolator if --protein T is set.", true);
  InitBoolParam("no-prune-proteins", false,
    "Peptides with low score will not be pruned before calculating protein probabilities "
    "(--protein T must be set).",
    "Available for crux percolator if --protein T is set.", true);
  InitIntParam("deepness", 0, 0, 2,
    "Setting deepness 0, 1, or 2 from low depth to high depth (less computational time) "
    "of the grid search for estimation Alpha,Beta and Gamma parameters for fido "
    "(--protein T must be set). Default value is 0.",
    "Available for crux percolator if --protein T is set.", true);
  InitBoolParam("original-output", false,
    "Output the standalone Percolator tab-delimited output.",
    "Available for crux percolator.", false);
  // **** Tide arguments ****
  InitArgParam("spectrum records file",
    "A spectrum records file generated by a previous run of crux tide-search "
    "using the store-spectra parameter.");
  InitArgParam("tide spectra file",
    "The name of the file from which to parse the fragmentation spectra, in any "
    "of the file formats supported by ProteoWizard. Alternatively, the argument "
    "may be a binary spectrum file produced by a previous run of crux "
    "tide-search using the store-spectra parameter.");
  InitArgParam("tide database index",
    "A directory containing a database index created by a previous run of crux "
    "tide-index.");
  // **** Tide options ****
  InitStringParam("decoy-format", "shuffle",
    "Include a decoy version of every peptide by shuffling or reversing the "
    "target sequence or protein. <string>=none|shuffle|peptide-reverse|"
    "protein-reverse. Default=shuffle.",
    "Available for tide-index", true);
  InitBoolParam("monoisotopic-precursor", true,
    "Use monoisotopic precursor masses rather than average mass for precursor. "
    "Default=T.",
    "Available for tide-index", true);
  InitStringParam("mods-spec", "C+57.02146",
    "Expression for static and variable mass modifications to include. "
    "Specify a comma-separated list of modification sequences of the form: "
    "C+57.02146,2M+15.9949,1STY+79.966331,...",
    "Available for tide-index", true);
  InitStringParam("cterm-peptide-mods-spec", "",
    "Specifies C-terminal static and variable mass modifications on peptides. "
    "Specify a comma-separated list of C-terminal modification sequences of the form: "
    "X+21.9819",
    "Available for tide-index", true);
  InitStringParam("nterm-peptide-mods-spec", "",
    "Specifies N-terminal static and variable mass modifications on peptides. "
    "Specify a comma-separated list of N-terminal modification sequences of the form: "
    "1E-18.0106,C-17.0265",
    "Available for tide-index", true);
  InitStringParam("cterm-protein-mods-spec", "",
    "Specifies C-terminal static and variable mass modifications on proteins. "
    "Specify a comma-separated list of C-terminal protein modification sequences of the form: "
    ",...",
    "Available for tide-index", true);
  InitStringParam("nterm-protein-mods-spec", "",
    "Specifies N-terminal static and variable mass modifications on proteins. "
    "Specify a comma-separated list of N-terminal protein modification sequences of the form: "
    ",...",
    "Available for tide-index", true);
  InitStringParam("store-spectra", "",
    "Specify the name of the file where the binarized fragmentation spectra "
    "will be stored.",
    "Available for tide-search", true);
  InitBoolParam("exact-p-value", false,
    "Uses exact P-value calculation for peptide-spectrum-matching. "
    "Default=F.",
    "Available for tide-search", true);
  InitBoolParam("concat", false,
    "Output target and decoy PSMs into a single file.",
    "Available for tide-search", true);
  InitBoolParam("file-column", true,
    "Include the file column in tab-delimited output.",
    "Available for tide-search", true);
  // Same as remove_precursor_peak and remove_precursor tolerance in Comet
  InitBoolParam("remove-precursor-peak", false,
    "Remove peaks around the precursor m/z.",
    "Available for tide-search.", true);
  InitDoubleParam("remove-precursor-tolerance", 1.5, 0, BILLION,
    "+- m/z tolerance for precursor peak removal. Default = 1.5.",
    "Available for print-processed spectra and tide-search.", true);
  InitBoolParam("clip-nterm-methionine", false,
    "This parameter controls whether Tide will automatically "
    "remove the N-terminal methionine from a sequence entry.",
    "Available for tide-index.", true);
  InitBoolParam("use-neutral-loss-peaks", false,
    "Controls whether neutral loss ions are considered in the search. "
    "Two types of neutral losses are included and are applied only to "
    "singly charged b- and y-ions: loss of ammonia (NH3, 17.0086343 Da) "
    "and H2O (18.0091422). Each neutral loss peak has intensity 1/5 of "
    "the primary peak",
    "Available for tide-search.", true);
  InitIntParam("max-precursor-charge", 5, 1, BILLION,
    "he maximum charge state of a spectra to consider in search. "
    "Default=5.",
    " Available for tide-search.", true);
  InitBoolParam("peptide-centric-search", false,
    "Carries out a peptide-centric search.  "
    "For each peptide the top N scoring spectra are reported,"
    "in contrary to the standard spectrum centric-search where top N "
    "scoring peptide are reported. Default = F.",
    "Available for tide-search.", true);
  InitIntParam("elution-window-size", 0, 0, 10,
    "Size of the elution window used in smoothing score in DIA mode. "
    "Used only with peptide-centric-search if greater than 0. A score of a psms "
    "centred in the window is substituted by the geometric mean of the scores "
    "in the window. If windows size is even, then it is increased by 1. Default = 0.",
    "Available for tide-search.", false);
  InitBoolParam("skip-decoys", true,
    "Skips decoys when reading a Tide index. Default = true.",
    "Available for read-tide-index", false);

  /*
   * Comet parameters
   */
  InitArgParam("input spectra",
    "The name of file (in MS2 format) from which to parse the spectra.");
  InitArgParam("database_name",
    "A full or relative path to the sequence database, "
    "in FASTA format, to search. Example databases include "
    "RefSeq or UniProt.  Database can contain amino acid "
    "sequences or nucleic acid sequences. If sequences are "
    "amino acid sequences, set the parameter \"nucleotide_reading_frame = 0\". "
    "If the sequences are nucleic acid sequences, you must instruct Comet to "
    "translate these to amino acid sequences. Do this by setting "
    "nucleotide_reading_frame\" to a value between 1 and 9.");
  InitIntParam("decoy_search", 0, 0, 2,
    "0=no (default), 1=concatenated search, 2=separate search",
    "option for Comet only", true);
  InitIntParam("num_threads", 0, 0, 32, 
    "0=poll CPU to set num threads; else specify num threads directly (max 32)",
    "option for Comet only", true);
  InitStringParam("output_suffix", "",
    "specifies the suffix string that is appended to the base output name "
    "for the pep.xml, pin.xml, txt and sqt output files."
    "Default = \"\"",
    "Available for comet.",true);
  InitDoubleParam("peptide_mass_tolerance", 3.0, 0, BILLION,
    "Controls the mass tolerance value.  The mass tolerance "
    "is set at +/- the specified number i.e. an entered value "
    "of \"1.0\" applies a -1.0 to +1.0 tolerance. "
    "The units of the mass tolerance is controlled by the parameter "
    "\"peptide_mass_units\". ", 
    "option for Comet only",true);
  InitIntParam("peptide_mass_units", 0, 0, 2,
    "0=amu, 1=mmu, 2=ppm",
    "option for Comet only", true);
  InitIntParam("mass_type_parent", 1, 0, 1,
    "0=average masses, 1=monoisotopic masses","option for Comet only", true);
  InitIntParam("mass_type_fragment", 1, 0, 1,
    "0=average masses, 1=monoisotopic masses","option for Comet only", true);
  InitIntParam("isotope_error", 0, 0, 2, 
    "0=off, 1=on -1/0/1/2/3 (standard C13 error), 2= -8/-4/0/4/8 (for +4/+8 labeling)",
    "option for Comet only", true);
  InitIntParam("search_enzyme_number", 1, 0, BILLION,
    "choose from list at end of this params file",
    "option for Comet only", true);
  InitIntParam("num_enzyme_termini", 2, 1, 2,
    "valid values are 1 (semi-digested), "
    "2 (fully digested, default), 8 N-term, 9 C-term",
    "option for Comet only", true);
  InitIntParam("allowed_missed_cleavage", 2, 0, 5,
    "maximum value is 5; for enzyme search",
    "option for Comet only", true);
  InitDoubleParam("fragment_bin_tol", 1.000507, 0, BILLION,
    "binning to use on fragment ions",
    "option for Comet only", true);
  InitDoubleParam("fragment_bin_offset", 0.40, 0, 1.0,
    "offset position to start the binning (0.0 to 1.0)",
    "option for Comet only", true);
  InitIntParam("theoretical_fragment_ions", 1, 0, 1,
    "0=default peak shape, 1=M peak only",
    "option for Comet only", true);
  InitIntParam("use_A_ions", 0, 0, 1, 
    "Controls whether or not A-ions are considered in the search (0 - no, 1 - yes)",
    "option for Comet only", true);
  InitIntParam("use_B_ions", 1, 0, 1, 
    "Controls whether or not B-ions are considered in the search (0 - no, 1 - yes)",
    "option for Comet only", true);
  InitIntParam("use_C_ions", 0, 0, 1, 
    "Controls whether or not C-ions are considered in the search (0 - no, 1 - yes)",
    "option for Comet only", true);
  InitIntParam("use_X_ions", 0, 0, 1, 
    "Controls whether or not X-ions are considered in the search (0 - no, 1 - yes)",
    "option for Comet only", true);
  InitIntParam("use_Y_ions", 1, 0, 1,
    "Controls whether or not Y-ions are considered in the search (0 - no, 1 - yes)",
    "option for Comet only", true);
  InitIntParam("use_Z_ions", 0, 0, 1, 
    "Controls whether or not Z-ions are considered in the search (0 - no, 1 - yes)",
    "option for Comet only", true);
  InitIntParam("use_NL_ions", 1, 0, 1,
    "0=no, 1= yes to consider NH3/H2O neutral loss peak",
    "option for Comet only", true);
  InitIntParam("use_sparse_matrix", 0, 0, 1,
    "Controls whether or not internal sparse matrix data representation is used.",
    "option for Comet only", true);
  InitIntParam("output_sqtfile", 0, 0, 1,
    "0=no, 1=yes  write sqt file",
    "option for Comet only", true);
  InitIntParam("output_pepxmlfile", 1, 0, 1,
    "0=no, 1=yes  write pep.xml file",
    "option for Comet only", true);
  InitIntParam("output_percolatorfile", 0, 0, 1,
    "0=no, 1=yes write percolator file",
     "option for Comet only", true);
  InitIntParam("output_txtfile", 1, 0, 1,
    "0=no, 1=yes  write tab-delimited text file",
    "option for Comet only (default 1)", true);
  InitIntParam("output_outfiles", 0, 0, 1,
    "0=no, 1=yes  write .out files",
    "option for Comet only", true);
  InitIntParam("print_expect_score", 1, 0, 1,
    "0=no, 1=yes to replace Sp with expect in out & sqt",
    "option for Comet.", true);
  InitIntParam("num_output_lines", 5, 1, BILLION,
    "num peptide results to show",
    "option for Comet.", true);
  InitIntParam("show_fragment_ions", 0, 0, 1,
    "0=no, 1=yes for out files only",
    "option for Comet.", true);
  InitIntParam("sample_enzyme_number", 1,0,10, 
    "Sample enzyme which is possibly different than the one applied to the search."
    "Used to calculate NTT & NMC in pepXML output (default=1 for trypsin).",
    "option for Comet. ", true);
  InitStringParam("scan_range", "0 0",
    "start and scan scan range to search; 0 as 1st entry "
    "ignores parameter",
    "option for Comet", true);
  InitStringParam("precursor_charge", "0 0",
    "precursor charge range to analyze; does not override "
    "mzXML charge; 0 as 1st entry ignores parameter",
    "option for Comet.", true);
  InitIntParam("ms_level", 2, 2, 3, 
    "MS level to analyze, valid are levels 2 (default) or 3",
    "option for Comet. ", true);
  InitStringParam("activation_method", "ALL", "ALL|CID|ECD|ETD|PQD|HCD|IRMPD",
    "<string>= ALL|CID|ECD|ETD|PQD|HCD|IRMPD. Default=All",
    "option for Comet. ", true);
  InitStringParam("digest_mass_range", "600.0 5000.0",
    "MH+ peptide mass range to analyze",
    "option for Comet.", true);
  InitIntParam("num_results", 50, 0, BILLION,
    "number of search hits to store internally",
    "option for Comet.", true);
  InitIntParam("skip_researching", 1, 0, 1,
    "for '.out' file output only, 0=search everything again "
    "(default), 1=don't search if .out exists",
    "option for Comet", true);
  InitIntParam("max_fragment_charge", 3, 1, 5,
    "set maximum fragment charge state to analyze (allowed max 5)",
    "option for Comet", true);
  InitIntParam("max_precursor_charge", 6, 1, 9,
    "set maximum precursor charge state to analyze (allowed max 9)",
    "option for Comet", true);
  InitIntParam("nucleotide_reading_frame", 0, 0, 9,
    "0=proteinDB, 1-6, 7=forward three, 8=reverse three, 9=all six",
    "option for Comet", true);
  InitIntParam("clip_nterm_methionine", 0, 0, 1,
    "0=leave sequences as-is; 1=also consider sequence w/o N-term methionine",
    "option for Comet", true);
  InitIntParam("spectrum_batch_size", 0, 0, BILLION,
    "max. # of spectra to search at a time; 0 to search the "
    "entire scan range in one loop",
    "option for Comet", true);
  InitIntParam("minimum_peaks", 10, 1, BILLION,
    "minimum num. of peaks in spectrum to search (default 10)",
    "option for Comet", true);
  InitDoubleParam("minimum_intensity", 0, 0, BILLION,
    "minimum intensity value to read in",
    "option for comet. ", true);
  InitIntParam("remove_precursor_peak", 0, 0, 2, 
    "0=no, 1=yes, 2=all charge reduced precursor peaks (for ETD)",
    "option for Comet. ", true);
  InitDoubleParam("remove_precursor_tolerance", 1.5, -BILLION, BILLION, 
    "+- Da tolerance for precursor removal",
    "option for Comet. ", true);
  InitStringParam("clear_mz_range", "0.0 0.0",
    "for iTRAQ/TMT type data; will clear out all peaks in the specified m/z range",
    "option for Comet", true);
  for (int i = 1; i <= 9; i++) {
    InitStringParam("variable_mod0" + StringUtils::ToString(i), "",
                    "Up to 9 variable modifications are supported\n"
                    "format: <mass> <residues> <0=variable/1=binary> <max mods per a peptide>\n"
                    "    e.g. 79.966331 STY 0 3",
                    "option for Comet", true);
  }
  InitIntParam("require_variable_mod", 0, 0, 1,
    "controls whether the analyzed peptides must contain at least one variable modification", 
    "option for Comet", true);
  InitIntParam("max_variable_mods_in_peptide", 5, 0, BILLION,
    "Specifies the total/maximum number of residues that can "
    "be modified in a peptide",
    "option for Comet", true);
  InitIntParam("override_charge", 0, 0, 1,
    "specifies the whether to override existing precursor charge state information when present "
    "in the files with the charge range specified by the \"precursor_charge\" parameter",
    "option for Comet", true);  
  InitDoubleParam("add_Cterm_peptide", 0, 0, BILLION,
    "Specifiy a static modification to the c-terminus of all peptides",
    "option for Comet", true);
  InitDoubleParam("add_Nterm_peptide", 0, 0, BILLION,
    "Specify a static modification to the n-terminus of all peptides",
    "option for Comet", true);
  InitDoubleParam("add_Cterm_protein", 0, 0, BILLION,
    "Specify a static modification to the c-terminal peptide of each protein",
    "option for Comet", true);
  InitDoubleParam("add_Nterm_protein", 0, 0, BILLION,
    "Specify a static modification to the n-terminal peptide of each protein",
    "option for Comet", true);
  for (char c = 'A'; c <= 'Z'; c++) {
    string aaString = string(1, c);
    string aaName = AminoAcidUtil::GetName(c);
    aaName = aaName.empty() ? "user_amino_acid" : StringUtils::Replace(aaName, " ", "_");
    InitDoubleParam("add_" + aaString + "_" + aaName,
                    c != 'C' ? 0 : CYSTEINE_DEFAULT, 0, BILLION,
                    "Specify a static modification to the residue " + aaString,
                    "option for Comet", true);
  }
  // **** q-ranker-barista arguments ****
  InitArgParam("database",
    "The program requires the FASTA format protein database files against "
    "which the search was performed. The protein database input may be a "
    "concatenated database or separate target and decoy databases; the "
    "latter is supported with the --separate-searches option, described "
    "below. In either case, Barista distinguishes between target and decoy "
    "proteins based on the presence of a decoy prefix on the sequence "
    "identifiers (see the --decoy-prefix option, below). The database can "
    "be provided in three different ways: (1) as a a single FASTA file "
    "with suffix \".fa\", \".fsa\" or \".fasta\", (2) as a text file "
    "containing a list of FASTA files, one per line, or (3) as a directory "
    "containing multiple FASTA files (identified via the filename suffixes "
    "\".fa\", \".fsa\" or \".fasta\")."); 
  InitArgParam("search results",
    "Q-ranker recognizes search results in tab-delimited format. Like the spectra, the "
    "search results can be provided as a single file, a list of files or a "
    "directory of files. Note, however, that the input mode for spectra and "
    "for search results must be the same; i.e., if you provide a list of "
    "files for the spectra, then you must also provide a list of files "
    "containing your search results. When the MS2 files and tab-delimited text files are "
    "provided via a file listing, Q-ranker assumes that the order of the MS2 "
    "files matches the order of the tab-delimited files. Alternatively, when the MS2 "
    "files and tab-delimited files are provided via directories, Q-ranker will search "
    "for pairs of files with the same root name but different extensions "
    "(\".ms2\" and \".txt\").");
  // **** q-ranker options. ****
  InitBoolParam("skip-cleanup", false, 
    "Q-ranker analysis begins with a pre-processsing step that creates a "
    "set of lookup tables which are then used during training. Normally, "
    "these lookup tables are deleted at the end of the Q-ranker analysis, "
    "but setting this option to T prevents the deletion of these tables. "
    "Subsequently, the Q-ranker analysis can be repeated more efficiently "
    "by specifying the --re-run option. Default = F.", 
    "Available for q-ranker and barista.", true);
  InitBoolParam("use-spec-features", true, 
    "Q-ranker uses enriched feature set derived from the spectra in ms2 "
    "files. It can be forced to use minimal feature set by setting the "
    "--use-spec-features option to F. Default T.", 
    "Available for q-ranker and barista.", true);
  InitStringParam("decoy_prefix", "decoy_",
    "Specifies the prefix of the protein names that indicates a decoy. "
    "Default = decoy_.",
    "Available for q-ranker and barista.", true);
  InitStringParam("re-run", "",
    "Re-run a previous Q-ranker analysis using a previously computed set of "
    "lookup tables.",
    "Available for q-ranker and barista.", true);
  InitStringParam("separate-searches", "",
    "If the target and decoy searches were run separately, rather than "
    "using a concatenated database, then Q-ranker will assume that the "
    "database search results provided as a required argument are from the "
    "target database search. This option then allows the user to specify "
    "the location of the decoy search results. Like the required arguments, "
    "these search results can be provided as a single file, a list of files "
    "or a directory. However, the choice (file, list or directory) must be "
    "consistent for the MS2 files and the target and decoy tab-delimited files. Also, "
    "if the MS2 and tab-delimited files are provided in directories, then Q-ranker "
    "will use the MS2 filename (foo.ms2) to identify corresponding target "
    "and decoy tab-delimited files with names like foo*.target.txt and "
    "foo*.decoy.txt. This naming convention allows the target and decoy txt "
    "files to reside in the same directory.",
    "Available for q-ranker and barista.", true);
  //**** Barista and QRanker options. ******
  InitBoolParam("list-of-files", false, 
    "Search result can be as a file or a list of files. This option"
    " allows users to specify the search results are provided as a list of files by " 
    "setting the --list-of-files option to T. Default= false.", 
    "Available for barista.",true);
  InitStringParam("optimization", "protein",
     "Specifies whether to do optimization at the protein, peptide or psm level. "
     "Default = protein.",
     "Available for barista.", true);
  /* analyze-matches parameter options */
  InitDoubleParam("pi-zero", 1.0, 0, 1, 
    "The estimated percent of target scores that are drawn from the "
    "null distribution.",
    "Used by assign-confidence, compute-q-values, percolator and q-ranker", false);
  InitBoolParam("smaller-is-better", false, 
    "Specify the semantics of the score, i.e., whether a smaller value "
    "implies a better match or vice versa.  Default is false. Specify this parameter "
    "T for \"exact p-value\" and F for \"xcorr score\". The default is F. ",
    "Used by assign-confidence.", true);
  InitStringParam("score", "xcorr score", 
    "Specify the column (for tab-delimited input) or tag (for XML input) "
    "used as input to the q-value estimation procedure. Default \"xcorr score\".",
    "Used by assign-confidence.", true);
  InitStringParam("estimation-method", "tdc", 
    "Specify the q-value calculation procedure to either "
    "target-decoy competition (tcd) or mix-max (mix-max). Default=tdc.",
    "Used by assign-confidence.", true);      
  InitStringParam("percolator-intraset-features", "F",
    "Set a feature for percolator that in later versions is not an option.",
    "Shouldn't be variable; hide from user.", false);
  // **** predict-peptide-ions options. ****
  InitStringParam("primary-ions", "by", "b|y|by|bya",
    "The ion series to predict (b,y,by,bya). Default='by' (both b and y ions).",
    "Only available for crux-predict-peptide-ions.  Set automatically to "
    "'by' for searching.", true);
  InitBoolParam("precursor-ions", false,
    "Predict the precursor ions, and all associated ions "
    "(neutral-losses, multiple charge states) consistent with the "
    "other specified options. (T,F) Default=F.",
    "Only available for crux-predict-peptide-ions.", true);
  InitIntParam("isotope", 0, 0, 2,
    "Predict the given number of isotope peaks (0|1|2). Default=0.",
    "Only available for crux-predict-peptide-ion.  Automatically set to "
    "0 for Sp scoring and 1 for xcorr scoring.", true);
  InitBoolParam("flanking", false, 
    "Predict flanking peaks for b and y ions (T,F). Default=F.",
    "Only available for crux-predict-peptide-ion.", true);
  InitStringParam("max-ion-charge", "peptide",
    "Predict theoretical ions up to max charge state (1,2,...,6) or up to the charge state "
    "of the peptide (peptide).  If the max-ion-charge is greater than the "
    "charge state of the peptide, then the max is the peptide charge. "
    "Default='peptide'.",
    "Available for predict-peptide-ions and search-for-xlinks. "
    "Set to 'peptide' for search.",
    true);
  InitIntParam("nh3", 0, -100, BILLION, 
    "Predict peaks with the given maximum number of nh3 neutral loss "
    "modifications. Default=0.",
    "Only available for crux-predict-peptide-ions.", true);
  InitIntParam("h2o", 0, -100, BILLION,
    "Predict peaks with the given maximum number of h2o neutral loss "
    "modifications. Default=0.",
    "Only available for crux-predict-peptide-ions.", true);
  // ***** spectral-counts aguments *****
  InitArgParam("input PSMs",
    "Name of file in text format which holds match results.");
  // also uses "protein-database"
  // ***** spectral-counts options *****
  InitStringParam("protein-database", "",
    "Fasta file of proteins.",
    "Option for spectral-counts", true);
  InitStringParam("input-ms2", "",
    "MS2 file corresponding to the psm file. Required for SIN.",
    "Available for spectral-counts with measure=SIN.", true);
  InitStringParam("threshold-type", "qvalue", "none|qvalue|custom",
    "What type of threshold to use when parsing matches "
    "none|qvalue|custom. Default=qvalue.",
    "used for crux spectral-counts", true);
  InitDoubleParam("threshold", 0.01,
    "The threshold to use for filtering matches. Default=0.01.",
    "Available for spectral-counts.  All PSMs with higher (or lower) than "
    "this will be ignored.", true);
  InitBoolParam("custom-threshold-min", true,
    "Direction of threshold for matches.  If true, then all matches "
    "whose value is <= threshold will be accepted.  If false, "
    "then all matches >= threshold will be accepted.  Used with "
    "custom-threshold parameter (default True).",
    "Available for spectral-counts.", true);
  InitStringParam("custom-threshold-name", "",
    "Use a name of custom threshold rather than (default NULL)",
    "Available for spectral-counts.", true);
  InitStringParam("measure", "NSAF", "RAW|NSAF|dNSAF|SIN|EMPAI",
    "Type of analysis to make on the match results: "
    "(RAW|NSAF|dNSAF|SIN|EMPAI). Default=NSAF. ", 
    "Available for spectral-counts.  RAW is raw counts, "
    "NSAF is Normalized Spectral Abundance Factor, "
    "dNSAF is Distributed Spectral Abundance Factor, "
    "SIN is Spectral Index Normalized and EMPAI is "
    "Exponentially Modified Protein Abundance Index", true);
  InitBoolParam("unique-mapping", false,
    "Ignore peptides with multiple mappings to proteins (T,F). Default=F.",
    "Available for spectral-counts.", true);
  InitStringParam("quant-level", "protein", "protein|peptide",
    "Quantification at protein or peptide level (protein,peptide). Default=protein.",
    "Available for spectral-counts and either NSAF and SIN.", true);
  InitStringParam("parsimony", "none", "none|simple|greedy",
    "Perform parsimony analysis on the proteins and report "
    "a parsimony rank column in output file. "
    "Default=none. Can be <string>=none|simple|greedy",
    "Available for spectral-counts.", true);
  InitBoolParam("mzid-use-pass-threshold", false,
    "Use mzid's passThreshold attribute to filter matches. Default false.",
    "Used when parsing mzIdentML files.", true);
  // ***** static mods *****
  for (char c = 'A'; c <= 'Z'; c++) {
    double deltaMass = (c != 'C') ? 0 : CYSTEINE_DEFAULT;
    bool visible = (c != 'B' && c != 'J' && c != 'O' && c != 'U' && c != 'X' && c != 'Z');
    InitDoubleParam(string(1, c), deltaMass,
      "Change the mass of all amino acids '" + string(1, c) + "' by the "
      "given amount.", "Default=" + StringUtils::ToString(deltaMass) + ".",
      visible);
  }
  /* get-ms2-spectrum options */
  InitBoolParam("stats", false, 
    "Print to stdout additional information about the spectrum.",
    "Avaliable only for crux-get-ms2-spectrum.  Does not affect contents "
    "of the output file.", true);
  // **** xlink-predict-peptide-ions options ****
  InitStringParam("peptide A", "", 
    "The sequence of peptide A.",
    "Argument for xlink-predict-peptide-ions.", false);
  InitStringParam("peptide B", "", 
    "The sequence of peptide B.",
    "Argument for xlink-predict-peptide-ions.", false);
  InitIntParam("pos A", 0 , 0, BILLION, 
    "Position of xlink on peptide A",
    "Available for xlink-predict-peptide-ions.", false);
  InitIntParam("pos B", 0 , 0, BILLION, 
    "Position of xlink on peptide B",
    "Available for xlink-predict-peptide-ions.", false);
  InitBoolParam("print-theoretical-spectrum", false,
    "Print the theoretical spectrum",
    "Available for xlink-predict-peptide-ions (Default=F).", true);
  InitBoolParam("use-old-xlink", true /* Turn to false later */,
    "Use old xlink searching algorithm",
    "Available for search-for-xlinks program (Default=T).", false);
  // **** xlink-score-spectrum options ****
  InitStringParam("xlink-score-method", "composite", "composite|modification|concatenated",
    "Score method for xlink {composite, modification, concatenated}. Default=composite.",
    "Argument for xlink-score-spectrum.", false);
  // **** search-xlink options ****
  InitStringParam("isotope-windows", "0",
    "List of integers of isotopic masses to search",
    "Used for crux search-for-xlinks", true);
  InitBoolParam("xlink-print-db", false,
    "Print the database in tab delimited format to xlink_peptides.txt",
    "Used for testing the candidate generatation (Default=F).", false);
  InitBoolParam("xlink-include-linears", true, 
    "Include linear peptides in the "
    "database.  Default=T.",
    "Available for crux search-for-xlinks program (Default=T).", true);
  InitBoolParam("xlink-include-deadends", true, 
    "Include dead-end peptides in the "
    "database.  Default=T.",
    "Available for crux search-for-xlinks program.", true);
  InitBoolParam("xlink-include-selfloops", true, 
    "Include self-loop peptides in the database.  Default=T.",
    "Available for crux search-for-xlinks program.", true);
  InitStringParam("xlink-prevents-cleavage", "K",
    "List of amino acids that xlinker can prevent cleavage",
    "Available for search-for-xlinks program (Default=K).",
    false /*TODO - turn this to true after new xlink code is implemented */);
  InitDoubleParam("precursor-window-weibull", 20.0, 0, 1e6, 
    "Search decoy-peptides within +/- "
    " 'mass-window-decoy' of the spectrum mass.  Default=20.0.",
    "Available for crux search-for-xlinks. ",
    true);
  InitStringParam("precursor-window-type-weibull", "mass", "mass|mz|ppm",
    "Window type to use for selecting "
    "decoy peptides from precursor mz. <string>=mass|mz|ppm. Default=mass.",
    "Available for crux search-for-xlinks", true);
  InitArgParam("link sites",
    "Comma delimited pair of amino acid link sites, e.g., A:K,A:D.");
  InitArgParam("link mass",
    "Mass modification of a cross link between two amino acids.");
  InitIntParam("min-weibull-points", 4000, 1, BILLION, 
    "Minimum number of points for estimating the "
    "Weibull parameters.  Default=4000.",
    "Available for crux search-for-xlinks", true);
  InitIntParam("max-xlink-mods", 0, 0, BILLION,
    "Maximum number of modifications allowed on a crosslinked peptide. Default=0.",
    "Available for crux search-for-xlinks", true);
  /* hardklor parameters */
  InitStringParam("hardklor-algorithm", "fast-fewest-peptides",
    "basic|fewest-peptides|fast-fewest-peptides|fewest-peptides-choice|"
    "fast-fewest-peptides-choice", 
    "Choose the algorithm for analyzing combinations of "
    "multiple peptide or protein isotope distributions. "
    "(basic | fewest-peptides | fast-fewest-peptides | "
    "fewest-peptides-choice | fast-fewest-peptides-choice) "
    "Default=fast-fewest-peptides.", "Available for crux hardklor", true);
  InitStringParam("cdm", "Q",
    "Choose the charge state determination method. (B|F|P|Q|S). Default=Q.",
    "Available for crux hardklor", true);
  InitIntParam("min-charge", 1, 1, BILLION,
    "Set the minimum charge state to look for when analyzing a spectrum. Default=1.",
    "Available for crux hardklor", true);
  InitIntParam("max-charge", 5, 1, BILLION,
    "Set the maximum charge state to look for when analyzing a spectrum. Default=5.",
    "Available for crux hardklor", true);
  InitDoubleParam("corr", 0.85, 0, 1.0, 
    "Set the correlation threshold [0,1.0] to accept a predicted "
    "isotope distribution.  Default=0.85",
    "Available for crux hardklor", true);
  InitIntParam("depth", 3, 1, BILLION,
    "Set the depth of combinatorial analysis. Default 3.",
    "Available for crux hardklor", true);
  InitBoolParam("distribution-area", false,
    "Reports peptide intensities as the distribution area. Default false.",
    "Available for crux hardklor", true);
  InitStringParam("averagine-mod", "",
    "Include alternative averagine models in the analysis that  "
    "incorporate additional atoms or isotopic enrichments.",
    "Available for crux hardklor", true);
  InitStringParam("mzxml-filter", "none",
    "Set a filter for mzXML files. Default=none",
    "Available for crux hardklor", true);
  InitBoolParam("no-base", false,
    "Specify \"no base\" averagine. Only modified averagine models "
    "will be used in the analysis. Default = F ",
    "Available for crux hardklor", true);
  InitIntParam("max-p", 10, 1, BILLION,
    "Set the maximum number of peptides or proteins that are "
    "estimated from the peaks found in a spectrum segment. The "
    "default value is 10.",
    "Available for crux hardklor", true);  
  InitDoubleParam("resolution", 100000, 1, BILLION,
    "Set the resolution of the observed spectra at m/z 400. "
    "Used in conjunction with --instrument The default is 100000.",
    "Available for crux hardklor", true);
  InitBoolParam("centroided", false,
    "Are spectra centroided?  Default false.",
    "Available for crux hardklor", true);
  InitStringParam("instrument", "fticr", "fticr|orbitrap|tof|qit",
    "Type of instrument (fticr|orbitrap|tof|qit) on which the data was "
    "collected. Used in conjuction with --resolution. The default is fticr.",
    "Available for crux hardklor", true);
  InitIntParam("sensitivity", 2, 0, 3,
    "Set the sensitivity level. There are four levels, 0 (low), 1 (moderate), "
    "2 (high), and 3 (max). The default value is 2.",
    "Available for crux hardklor", true);
  InitDoubleParam("signal-to-noise", 1.0, 0.0, BILLION,
    "Set the signal-to-noise threshold. Any integer or decimal "
    "value greater than or equal to 0.0 is valid. The default value is 1.0.",
    "Available for crux hardklor", true);
  InitDoubleParam("sn-window", 250.0, 0.0, BILLION,
    "Set the signal-to-noise window length (in m/z). Because noise may "
    "be non-uniform across a spectra, this value adjusts the segment size "
    "considered when calculating a signal-over-noise ratio. The default "
    "value is 250.0.",
    "Available for crux hardklor", true);
  InitBoolParam("static-sn", true,
    "If true, Hardklor will calculate the local noise levels across the "
    "spectrum using --sn-window, then select a floor of this set of noise "
    "levels to apply to the whole spectrum.",
    "Available for crux hardklor", true);
  InitStringParam("mz-window", "",
    "Restrict analysis to only a small window in each segment ( (min-max) in m/z). "
    "The user must specify the starting and ending m/z values between which "
    "the analysis will be performed. By default the whole spectrum is analyzed.",
    "Available for crux hardklor", true);
  InitDoubleParam("max-width", 4.0, 0.0, BILLION,
    "Set the maximum width of any set of peaks in a spectrum when computing the "
    "results (in m/z). Thus, if the value was 5.0, then sets of peaks greater "
    "than 5 m/z are divided into smaller sets prior to analysis. The default value is 4.0.",
    "Available for crux hardklor", true);
  InitStringParam("hardklor-options", "", 
    "Directly set hardklor options",
    "Available for crux hardklor", false);
  /* bullseye parameters */
  InitArgParam("MS1 spectra",
    "The name of a file from which to parse high-resolution spectra of intact peptides. "
    "The file may be in MS1 (.ms1), binary MS1 (.bms1), compressed MS1 (.cms1), or "
    "mzXML (.mzXML) format. ");
  InitArgParam("MS2 spectra",
    "The name of a file from which to parse peptide fragmentation spectra. The file may "
    "be in MS2 (.ms2), binary MS2 (.bms2), compressed MS2 (.cms2) or mzXML (.mzXML) format.");
  InitStringParam("hardklor-file", "",
    "Input hardklor file into bullseye",
    "Hidden option for crux bullseye.", false);
  InitDoubleParam("max-persist", 2.0, 0, BILLION,
    "Ignore peptides that persist for this length. The unit of time is whatever unit is "
    "used in your data file (usually minutes). These peptides are considered contaminants." 
    " Default = 2.0.",
    "Available for crux bullseye", true);
  InitBoolParam("exact-match", false, 
    "Require an exact match to the precursor ion. Rather than use wide precursor boundaries, "
    "this flag forces Bullseye to match precursors to the base isotope peak identified in "
    "Hardklor. The tolerance is set with the --persist-tolerance flag. Default = F.",
    "Available for crux bullseye", true);
  InitIntParam("gap-tolerance", 1, 0, BILLION,
    "Gap size tolerance when checking for peptides across consecutive MS1 scans. Used in "
    "conjunction with --scan-tolerance. Default = 1.",
    "Available for crux bullseye", true);
  InitDoubleParam("bullseye-min-mass", 600, 0, BILLION,
    "The minimum mass of peptides to consider. Default=600.",
    "Available from command line or parameter file for crux bullseye", true);
  InitDoubleParam("bullseye-max-mass", 8000, 1, BILLION, 
    "The maximum mass of peptides to consider. Default=8000.",
    "Available from command line or parameter file for crux bullseye", true);
  InitDoubleParam("exact-tolerance", 10.0, 0, BILLION,
    "Set the tolerance (+/-ppm) for exact match searches. Default = 10.0.",
    "Available for crux bullseye", true);
  InitDoubleParam("persist-tolerance", 10.0, 0, BILLION,
    "Set the tolerance (+/-ppm) for finding persistent peptides. Default = 10.0.",
    "Available for crux bullseye", true);
  InitIntParam("scan-tolerance", 3, 0, BILLION,
    "Number of consecutive MS1 scans over which a peptide must be observed to "
    "be considered real. Gaps in persistence are allowed when setting --gap-tolerance. "
    "Default = 3.",
    "Available for crux bullseye", true);
  InitDoubleParam("retention-tolerance", 0.5, 0, BILLION,
    "Set the tolerance (+/-units) around the retention time over which a peptide "
    "can be matched to the MS/MS spectrum. The unit of time is whatever unit is "
    "used in your data file (usually minutes). Default = 0.5.",
    "Available for crux bullseye", true);
  InitStringParam("spectrum-format", "",
    "The format to write the output spectra to. By default, the spectra will be "
    "output in the same format as the MS/MS input.",
    "Available for crux bullseye", true);
  /* crux-util parameters */
  InitBoolParam("ascending", true,
    "Sort in ascending order.  Otherwise, descending. "
    "Default: True.",
    "Available for sort-by-column", true);
  InitArgParam("tsv file",
    "Path to a delimited file (-) for standard input");
  InitStringParam("delimiter", "tab",
    "Character delimiter to use when parsing a delimited file (Default tab).",
    "Available for the delimited utility programs.", false);
  InitArgParam("column names",
    "List of column names separated by a comma");
  InitArgParam("column name",
    "Name of the column to do the operation on");
  InitArgParam("column value",
    "value of the column");
  InitBoolParam("header", true,
    "Print the header line of the tsv file. Default=T.",
    "Available for crux extract-columns and extract-rows", true);
  InitStringParam("column-type", "string", "int|real|string",
    "Specifies the data type the column contains "
    "(int|real|string) Default: string",
    "Available for crux extract-rows", true);
  InitStringParam("comparison", "eq", "eq|gt|gte|lt|lte|neq",
    "Specifies the operator that is used to compare an "
    "entry in the specified column to the value given "
    "on the command line.  (eq|gt|gte|lt|lte|neq). "
    "Default: eq.",
    "Available for crux extract-rows", true);
  // create-docs
  InitArgParam("tool name",
    "Specifies the tool to generate documentation for. If value is 'list', "
    "a list of available tools will be given.");
  InitStringParam("doc-template", "",
    "Specifies the main template to be used for generating documentation.",
    "Available for crux create-docs", false);
  InitStringParam("doc-input-template", "",
    "Specifies the template to be used for inputs when generating "
    "documentation.",
    "Available for crux create-docs", false);
  InitStringParam("doc-output-template", "",
    "Specifies the template to be used for outputs when generating "
    "documentation.",
    "Available for crux create-docs", false);
  InitStringParam("doc-option-category-template", "",
    "Specifies the template to be used for option categories when generating "
    "documentation.",
    "Available for crux create-docs", false);
  InitStringParam("doc-option-template", "",
    "Specifies the template to be used for options when generating "
    "documentation.",
    "Available for crux create-docs", false);
}

void Params::Categorize() {
  set<string> items;

  items.clear();
  items.insert("max-persist");
  items.insert("persist-tolerance");
  items.insert("scan-tolerance");
  items.insert("gap-tolerance");
  items.insert("bullseye-max-mass");
  items.insert("bullseye-min-mass");
  container_.AddCategory("Identifying PPIDs in MS1 spectra", items);

  items.clear();
  items.insert("exact-match");
  items.insert("exact-tolerance");
  items.insert("retention-tolerance");
  container_.AddCategory("Matching PPIDs to MS2 spectra", items);

  items.clear();
  items.insert("max-length");
  items.insert("min-length");
  items.insert("max-mass");
  items.insert("min-mass");
  items.insert("monoisotopic-precursor");
  items.insert("clip-nterm-methionine");
  container_.AddCategory("Peptide properties", items);

  items.clear();
  items.insert("mods-spec");
  items.insert("nterm-peptide-mods-spec");
  items.insert("cterm-peptide-mods-spec");
  items.insert("max-mods");
  items.insert("mod");
  for (char c = 'A'; c <= 'Z'; c++) {
    items.insert(string(1, c));
  }
  container_.AddCategory("Amino acid modifications", items);

  items.clear();
  items.insert("decoy-format");
  items.insert("keep-terminal-aminos");
  container_.AddCategory("Decoy database generation", items);

  items.clear();
  items.insert("enzyme");
  items.insert("custom-enzyme");
  items.insert("digestion");
  items.insert("missed-cleavages");
  container_.AddCategory("Enzymatic digestion", items);

  items.clear();
  items.insert("max-precursor-charge");
  items.insert("max-ion-charge");
  items.insert("peptide-centric-search");
  items.insert("exact-p-value");
  items.insert("precursor-window");
  items.insert("precursor-window-type");
  items.insert("compute-sp");
  items.insert("spectrum-min-mz");
  items.insert("spectrum-max-mz");
  items.insert("min-peaks");
  items.insert("spectrum-charge");
  items.insert("scan-number");
  items.insert("remove-precursor-peak");
  items.insert("remove-precursor-tolerance");
  items.insert("print-search-progress");
  items.insert("use-flanking-peaks");
  items.insert("use-neutral-loss-peaks");
  items.insert("mz-bin-width");
  items.insert("mz-bin-offset");
  items.insert("precursor-window-weibull");
  items.insert("precursor-window-type-weibull");
  items.insert("mod-mass-format");
  items.insert("fragment-mass");
  items.insert("isotopic-mass");
  items.insert("isotope-windows");
  items.insert("compute-p-values");
  container_.AddCategory("Search parameters", items);

  items.clear();
  items.insert("protein");
  items.insert("alpha");
  items.insert("beta");
  items.insert("gamma");
  items.insert("allow-protein-group");
  items.insert("protein-level-pi0");
  items.insert("empirical-protein-q");
  items.insert("group-proteins");
  items.insert("no-prune-proteins");
  items.insert("deepness");
  container_.AddCategory("Fido options", items);

  items.clear();
  items.insert("spectrum-format");
  items.insert("spectrum-parser");
  items.insert("top-match");
  items.insert("store-spectra");
  items.insert("xlink-print-db");
  items.insert("fileroot");
  items.insert("output-dir");
  items.insert("overwrite");
  items.insert("txt-output");
  items.insert("sqt-output");
  items.insert("pepxml-output");
  items.insert("mzid-output");
  items.insert("pin-output");
  items.insert("pout-output");
  items.insert("parameter-file");
  items.insert("verbosity");
  container_.AddCategory("Input and output", items);
}

bool Params::GetBool(const string& name) {
  return GetParam(name)->GetBool();
}

int Params::GetInt(const string& name) {
  return GetParam(name)->GetInt();
}

double Params::GetDouble(const string& name) {
  return GetParam(name)->GetDouble();
}

string Params::GetString(const string& name) {
  return GetParam(name)->GetString();
}

const vector<string>& Params::GetStrings(const string& name) {
  Param* param = GetParam(name);
  if (!param->IsArgument()) {
    throw runtime_error("Parameter '" + name + "' is not an argument");
  }
  return ((ArgParam*)param)->GetStrings();
}

string Params::GetUsage(const string& name) {
  return GetParam(name)->GetUsage();
}

string Params::GetFileNotes(const string& name) {
  return GetParam(name)->GetFileNotes();
}

bool Params::IsVisible(const string& name) {
  return GetParam(name)->IsVisible();
}

bool Params::IsArgument(const string& name) {
  return GetParam(name)->IsArgument();
}

string Params::GetType(const string& name) {
  return GetParam(name)->GetType();
}

bool Params::IsDefault(const string& name) {
  return GetParam(name)->IsDefault();
}

bool Params::Exists(const string& name) {
  return container_.Get(name) != NULL;
}

void Params::Set(const string& name, bool value) {
  container_.CanModifyCheck();
  Param* param = GetParam(name);
  param->Set(value);
  param->ThrowIfInvalid();
}

void Params::Set(const string& name, int value) {
  container_.CanModifyCheck();
  Param* param = GetParam(name);
  param->Set(value);
  param->ThrowIfInvalid();
}

void Params::Set(const string& name, double value) {
  container_.CanModifyCheck();
  Param* param = GetParam(name);
  param->Set(value);
  param->ThrowIfInvalid();
}

void Params::Set(const string& name, const char* value) {
  Set(name, string(value));
}

void Params::Set(const string& name, const string& value) {
  container_.CanModifyCheck();
  Param* param = GetParam(name);
  param->Set(value);
  param->ThrowIfInvalid();
}

void Params::AddArgValue(const string& name, const string& value) {
  container_.CanModifyCheck();
  Param* param = GetParam(name);
  if (!param->IsArgument()) {
    throw runtime_error("Cannot add value to '" + name + "', it is not an argument");
  }
  ((ArgParam*)param)->AddValue(value);
}

void Params::Finalize() {
  if (container_.Finalized()) {
    return;
  }

  if (GetString("enzyme") == "no-enzyme") {
    Set("digestion", "non-specific-digest");
    Set("missed-cleavages", 500);
  }

  for (char c = 'A'; c <= 'Z'; c++) {
    double deltaMass = GetDouble(string(1, c));
    increase_amino_acid_mass(c, deltaMass);
  }

  translate_decoy_options();

  string customEnzyme = GetString("custom-enzyme");
  if (!customEnzyme.empty()) {
    parse_custom_enzyme(customEnzyme);
    Set("enzyme", "custom-enzyme");
  }

  if (GetString("enzyme") == "no-enzyme") {
    Set("digestion", "non-specific-digest");
  } else if (GetString("digestion") == "non-specific-digest") {
    Set("enzyme", "no-enzyme");
  }

  double new_value = GetDouble("mz-bin-width");
// ***************************
#ifdef _MSC_VER
  // Peculiarities of Windows floating point handling 
  // results in us getting 0.0 here rather than Nan
  // FIXME: is there a more portable way of checking
  // that a floating point value has not been set?
  if (new_value == 0.0) {
#else
  if (isnan(new_value)) {
#endif
    // If no width specified, choose based on mass type.
    if (get_mass_type_parameter("fragment-mass") == MONO) {
      new_value = BIN_WIDTH_MONO;
    } else {
      new_value = BIN_WIDTH_AVERAGE;
    }

    Set("mz-bin-width", new_value);
  }
// ***************************

  container_.Finalize();
}

void Params::Write(ofstream* file) {
  if (file == NULL || !file->good()) {
    throw runtime_error("Bad file stream for writing parameter file");
  }

  *file << "# comet_version 2015.01 rev. 0" << endl
        << "# Comet MS/MS search engine parameters file." << endl
        << "# Everything following the \'#\' symbol is treated as a comment." << endl
        << endl;

  for (vector<const Param*>::const_iterator i = Begin(); i != End(); i++) {
    string name = (*i)->GetName();
    // Print mods and Comet parameters later
    if (!(*i)->IsVisible() ||
        name == "mod" || name == "cmod" || name == "nmod" ||
        name.find('_') != string::npos) {
      continue;
    }
    *file << (*i)->GetParamFileString() << endl;
  }

  print_mods_parameter_file(file, "mod", get_aa_mod_list);
  print_mods_parameter_file(file, "nmod", get_n_mod_list);
  print_mods_parameter_file(file, "cmod", get_c_mod_list);

  // Print Comet parameters
  *file << "#################" << endl
        << "#Comet Parameters" << endl
        << "#################" << endl;
  for (vector<const Param*>::const_iterator i = Begin(); i != End(); i++) {
    string name = (*i)->GetName();
    // Print mods and Comet parameters later
    if (!(*i)->IsVisible() ||
        name == "mod" || name == "cmod" || name == "nmod" ||
        name.find('_') == string::npos) {
      continue;
    }
    *file << (*i)->GetParamFileString() << endl;
  }

  *file << "#" << endl
        << "# COMET_ENZYME_INFO _must_ be at the end of this parameters file" << endl
        << "#" << endl
        << "[COMET_ENZYME_INFO]" << endl;

  const vector<string>& cometEnzymes = get_comet_enzyme_info_lines();
  if (cometEnzymes.empty()) {
    *file << "0.  No_enzyme                      0  -          -" << endl
          << "1.  Trypsin                        1  KR         P" << endl
          << "2.  Trypsin/P                      1  KR         -" << endl
          << "3.  Lys_C                          1  K          P" << endl
          << "4.  Lys_N                          0  K          -" << endl
          << "5.  Arg_C                          1  R          P" << endl
          << "6.  Asp_N                          0  D          -" << endl
          << "7.  CNBr                           1  M          -" << endl
          << "8.  Glu_C                          1  DE         P" << endl
          << "9.  PepsinA                        1  FL         P" << endl
          << "10. Chymotrypsin                   1  FWYL       P" << endl;
    /*TODO: Put these back in after we figure out what to do with enzyme info
    *file << "11. Elastase                       1  ALIV       P" << endl
          << "12. Clostripain                    1  R          -" << endl
          << "13. Iodosobenzoate                 1  W          -" << endl
          << "14. Proline_Endopeptidase          1  P          -" << endl
          << "15. Staph_Protease                 1  E          -" << endl
          << "16. Modified_Chymotrypsin          1  FWYL       P" << endl
          << "17. Elastase_Trypsin_Chymotrypsin  1  ALIVKRWFY  P" << endl;
    */
  } else {
    for (vector<string>::const_iterator i = cometEnzymes.begin();
         i != cometEnzymes.end();
         i++) {
      *file << *i << endl;
    }
  }
}

map<string, Param*>::const_iterator Params::BeginAll() {
  return container_.BeginAll();
}

map<string, Param*>::const_iterator Params::EndAll() {
  return container_.EndAll();
}

vector<const Param*>::const_iterator Params::Begin() {
  return container_.Begin();
}

vector<const Param*>::const_iterator Params::End() {
  return container_.End();
}

string Params::ProcessHtmlDocTags(string s, bool html) {
  // If html is true, instances of [[html:{text}]] become {text} and
  // instances of [[nohtml:{text}]] are removed.
  // If html is false, instances of [[nohtml:{text}]] become {text} and
  // instances of [[html:{text}]] are removed.

  const string OPEN_TAG = "[[";
  const string CLOSE_TAG = "]]";
  const string HTML_PREFIX = "html:";
  const string NO_HTML_PREFIX = "nohtml:";

  size_t pos, endPos = 0;
  while ((pos = s.find(OPEN_TAG, endPos)) != string::npos) {
    size_t prefixStart = pos + OPEN_TAG.length();
    if ((endPos = s.find(CLOSE_TAG, prefixStart)) == string::npos) {
      return s;
    }

    string fullOpen = OPEN_TAG;
    bool fullRemove;
    if (s.length() >= prefixStart + HTML_PREFIX.length() &&
        s.compare(prefixStart, HTML_PREFIX.length(), HTML_PREFIX) == 0) {
      fullOpen += HTML_PREFIX;
      fullRemove = !html;
    } else if (s.length() >= prefixStart + NO_HTML_PREFIX.length() &&
               s.compare(prefixStart, NO_HTML_PREFIX.length(), NO_HTML_PREFIX) == 0) {
      fullOpen += NO_HTML_PREFIX;
      fullRemove = html;
    } else {
      endPos = prefixStart;
      continue;
    }

    if (!fullRemove) {
      s.erase(pos, fullOpen.length());
      endPos -= fullOpen.length();
      s.erase(endPos, CLOSE_TAG.length());
    } else {
      s.erase(pos, endPos + CLOSE_TAG.length() - pos);
      endPos = pos;
    }
  }
  return s;
}

vector< pair< string, vector<string> > > Params::GroupByCategory(const vector<string>& options) {
  vector< pair< string, vector<string> > > groups;

  pair< string, vector<string> > uncategorizedPair = make_pair("", vector<string>(options));
  vector<string>& uncategorized = uncategorizedPair.second;

  const vector<ParamCategory>& categories = container_.GetCategories();
  // Iterate over all categories
  for (vector<ParamCategory>::const_iterator i = categories.begin();
       i != categories.end();
       i++) {
    bool any = false;
    // Iterate over each given option and check if it is in the category
    for (vector<string>::const_iterator j = options.begin(); j != options.end(); j++) {
      Param* param = GetParam(*j);
      // This option was in the category
      if (i->Items.find(param) != i->Items.end()) {
        if (!any) {
          any = true;
          groups.push_back(make_pair(i->Name, vector<string>()));
        }
        groups.back().second.push_back(*j);
        vector<string>::iterator iter;
        while ((iter = find(uncategorized.begin(), uncategorized.end(), *j)) !=
               uncategorized.end()) {
          uncategorized.erase(iter);
        }
      }
    }
  }

  if (!uncategorized.empty()) {
    groups.insert(groups.begin(), uncategorizedPair);
  }
  return groups;
}

void Params::InitBoolParam(
  const string& name,
  bool value,
  const string& usage,
  const string& fileNotes,
  bool visible
) {
  container_.Add(new BoolParam(name, usage, fileNotes, visible, value));
}

void Params::InitIntParam(
  const string& name,
  int value,
  int min,
  int max,
  const string& usage,
  const string& fileNotes,
  bool visible
) {
  container_.Add(new IntParam(name, usage, fileNotes, visible, value, min, max));
}

void Params::InitIntParam(
  const string& name,
  int value,
  const string& usage,
  const string& fileNotes,
  bool visible
) {
  container_.Add(new IntParam(name, usage, fileNotes, visible, value));
}

void Params::InitDoubleParam(
  const string& name,
  double value,
  double min,
  double max,
  const string& usage,
  const string& fileNotes,
  bool visible
) {
  container_.Add(new DoubleParam(name, usage, fileNotes, visible, value, min, max));
}

void Params::InitDoubleParam(
  const string& name,
  double value,
  const string& usage,
  const string& fileNotes,
  bool visible
) {
  container_.Add(new DoubleParam(name, usage, fileNotes, visible, value));
}

void Params::InitStringParam(
  const string& name,
  const string& value,
  const string& validValues,
  const string& usage,
  const string& fileNotes,
  bool visible
) {
  container_.Add(new StringParam(name, usage, fileNotes, visible, value,
                 StringUtils::Split(validValues, '|')));
}

void Params::InitStringParam(
  const string& name,
  const string& value,
  const string& usage,
  const string& fileNotes,
  bool visible
) {
  container_.Add(new StringParam(name, usage, fileNotes, visible, value));
}

void Params::InitArgParam(
  const string& name,
  const string& usage
) {
  container_.Add(new ArgParam(name, usage));
}

Param* Params::GetParam(const string& name) {
  Param* param = container_.Get(name);
  if (param == NULL) {
    throw runtime_error("Parameter '" + name + "' does not exist");
  }
  return param;
}

Params::Params() {
}

Params::~Params() {
}

// ***** Parameter container ***** //
Params::ParamContainer::ParamContainer()
  : finalized_(false) {
}

Params::ParamContainer::~ParamContainer() {
  for (map<string, Param*>::iterator i = params_.begin(); i != params_.end(); i++) {
     delete i->second;
  }
  //for (int i = 0; i < MAX_AA_MODS; i++) {
    //free_aa_mod(list_of_mods[i]);
    //list_of_mods[i] = NULL;
  //}
}

void Params::ParamContainer::Add(Param* param) {
  CanModifyCheck();
  param->ThrowIfInvalid();

  string paramName = param->GetName();
  if (!params_.insert(make_pair(paramName, param)).second) {
    throw runtime_error("Parameter '" + paramName + "' already exists");
  }
  if (!param->IsArgument()) {
    paramsOrdered_.push_back(param);
  }
}

Param* Params::ParamContainer::Get(const string& name) {
  map<string, Param*>::iterator i = params_.find(name);
  return (i == params_.end()) ? NULL : i->second;
}

bool Params::ParamContainer::Empty() const {
  return params_.empty();
}

bool Params::ParamContainer::Finalized() const {
  return finalized_;
}

map<string, Param*>::const_iterator Params::ParamContainer::BeginAll() const {
  return params_.begin();
}

map<string, Param*>::const_iterator Params::ParamContainer::EndAll() const {
  return params_.end();
}

vector<const Param*>::const_iterator Params::ParamContainer::Begin() const {
  return paramsOrdered_.begin();
}

vector<const Param*>::const_iterator Params::ParamContainer::End() const {
  return paramsOrdered_.end();
}

void Params::ParamContainer::Finalize() {
  finalized_ = true;
}

void Params::ParamContainer::CanModifyCheck() const {
  if (finalized_) {
    throw runtime_error("Parameters have been finalized and cannot be modified");
  }
}

void Params::ParamContainer::AddCategory(const string& name, const set<string>& params) {
  // Validate passed in set
  for (set<string>::const_iterator i = params.begin(); i != params.end(); i++) {
    if (Get(*i) == NULL) {
      throw runtime_error("Parameter '" + *i + "' does not exist");
    }
  }

  ParamCategory* category = NULL;
  // Check if this category already exists
  for (vector<ParamCategory>::iterator i = categories_.begin(); i != categories_.end(); i++) {
    if (i->Name == name) {
      category = &*i;
      break;
    }
  }

  // Create new category
  if (!category) {
    categories_.push_back(ParamCategory(name));
    category = &(categories_.back());
  }

  // Loop over parameters and add them to the category if they are in the passed in set
  for (vector<const Param*>::const_iterator i = Begin(); i != End(); i++) {
    string paramName = (*i)->GetName();
    if (params.find(paramName) != params.end()) {
      // Check if this parameter has already been categorized
      for (vector<ParamCategory>::const_iterator j = categories_.begin();
           j != categories_.end();
           j++) {
        if (j->Items.find(*i) != j->Items.end()) {
          throw runtime_error("Parameter '" + paramName + "' has already been categorized");
        }
      }
      // Add parameter to category
      category->Items.insert(*i);
    }
  }
}

const vector<Params::ParamCategory>& Params::ParamContainer::GetCategories() const {
  return categories_;
}

// ***** Parameter classes ***** //
//
// Param (base class)
//
Param::Param(const string& name,
             const string& usage,
             const string& fileNotes,
             bool visible)
  : name_(name), usage_(usage), fileNotes_(fileNotes), visible_(visible) {}
Param::~Param() {}
string Param::GetName() const { return name_; }
string Param::GetUsage() const { return usage_; }
string Param::GetFileNotes() const { return fileNotes_; }
bool Param::IsVisible() const { return visible_; }
bool Param::IsArgument() const { return false; }
void Param::ThrowIfInvalid() const {}
string Param::GetParamFileString() const {
  vector<string> lines =
    StringUtils::Split(Params::ProcessHtmlDocTags(usage_), '\n');
  vector<string> noteLines =
    StringUtils::Split(Params::ProcessHtmlDocTags(fileNotes_), '\n');
  lines.insert(lines.end(), noteLines.begin(), noteLines.end());
  stringstream ss;
  for (vector<string>::const_iterator i = lines.begin(); i != lines.end(); i++) {
    vector<string> formatted = StringUtils::Split(StringUtils::LineFormat(*i, 78), '\n');
    for (vector<string>::const_iterator j = formatted.begin(); j != formatted.end(); j++) {
      ss << "# " << *j << endl;
    }
  }
  ss << name_ << '=' << GetString() << endl;
  return ss.str();
}
void Param::Set(bool value) {
  throw runtime_error("Cannot set value of '" + name_ + "' from bool");
}
void Param::Set(int value) {
  throw runtime_error("Cannot set value of '" + name_ + "' from int");
}
void Param::Set(double value) {
  throw runtime_error("Cannot set value of '" + name_ + "' from double");
}
void Param::Set(const char* value) {
  Set(string(value));
}
void Param::Set(const string& value) {
  throw runtime_error("Cannot set value of '" + name_ + "' from string");
}
//
// BoolParam
//
BoolParam::BoolParam(const string& name,
                     const string& usage,
                     const string& fileNotes,
                     bool visible,
                     bool value)
  : Param(name, usage, fileNotes, visible), value_(value), original_(value) {}
string BoolParam::GetType() const { return "boolean"; }
bool BoolParam::IsDefault() const { return value_ == original_; }
bool BoolParam::GetBool() const { return value_; }
int BoolParam::GetInt() const { return value_ ? 1 : 0; }
double BoolParam::GetDouble() const { return value_ ? 1 : 0; }
string BoolParam::GetString() const { return value_ ? "true" : "false"; }
void BoolParam::Set(bool value) { value_ = value; }
void BoolParam::Set(int value) { value_ = value != 0; }
void BoolParam::Set(double value) { value_ = value != 0; }
void BoolParam::Set(const string& value) {
  try {
    value_ = FromString(value);
  } catch (...) {
    throw runtime_error("Invalid value for '" + name_ + "': " + "'" + value + "' "
                        "(expected boolean)");
  }
}
bool BoolParam::FromString(string s) {
  s = StringUtils::ToLower(s);
  if (s == "t" || s == "true") {
    return true;
  } else if (s == "f" || s == "false") {
    return false;
  }
  throw runtime_error("Cannot convert '" + s + "' to boolean");
}
//
// IntParam
//
IntParam::IntParam(const string& name,
                   const string& usage,
                   const string& fileNotes,
                   bool visible,
                   int value,
                   int min,
                   int max)
  : Param(name, usage, fileNotes, visible),
    value_(value), original_(value), min_(min), max_(max) {}
void IntParam::ThrowIfInvalid() const {
  if (value_ < min_ || value_ > max_) {
    throw runtime_error("Value of '" + name_ + "' must be between " +
                        StringUtils::ToString(min_) + " and " +
                        StringUtils::ToString(max_));
  }
}
string IntParam::GetType() const { return "integer"; }
bool IntParam::IsDefault() const { return value_ == original_; }
bool IntParam::GetBool() const { return value_ != 0; }
int IntParam::GetInt() const { return value_; }
double IntParam::GetDouble() const { return (double)value_; }
string IntParam::GetString() const { return StringUtils::ToString(value_); }
void IntParam::Set(bool value) { value_ = value ? 1 : 0; }
void IntParam::Set(int value) { value_ = value; }
void IntParam::Set(double value) { value_ = (int)value; }
void IntParam::Set(const string& value) {
  try {
    value_ = StringUtils::FromString<int>(value);
  } catch (...) {
    throw runtime_error("Invalid value for '" + name_ + "': " + "'" + value + "' "
                        "(expected int)");
  }
}
//
// DoubleParam
//
DoubleParam::DoubleParam(const string& name,
                         const string& usage,
                         const string& fileNotes,
                         bool visible,
                         double value,
                         double min,
                         double max)
  : Param(name, usage, fileNotes, visible),
    value_((FLOAT_T)value), original_((FLOAT_T)value), min_(min), max_(max) {}
    // TODO : For compatibility with old tests, convert to float
void DoubleParam::ThrowIfInvalid() const {
  if (value_ < min_ || value_ > max_) {
    throw runtime_error("Value of '" + name_ + "' must be between " +
                        StringUtils::ToString(min_) + " and " +
                        StringUtils::ToString(max_));
  }
}
string DoubleParam::GetType() const { return "float"; }
bool DoubleParam::IsDefault() const { return value_ == original_; }
bool DoubleParam::GetBool() const { return value_ != 0; }
int DoubleParam::GetInt() const { return (int)value_; }
double DoubleParam::GetDouble() const { return value_; }
string DoubleParam::GetString() const { return StringUtils::ToString(value_, 6); }
void DoubleParam::Set(bool value) { value_ = value ? 1 : 0; }
void DoubleParam::Set(int value) { value_ = (double)value; }
void DoubleParam::Set(double value) { value_ = value; }
void DoubleParam::Set(const string& value) {
  try {
    value_ = StringUtils::FromString<double>(value);
  } catch (...) {
    throw runtime_error("Invalid value for '" + name_ + "': " + "'" + value + "' "
                        "(expected float)");
  }
}
//
// StringParam
//
StringParam::StringParam(const string& name,
                         const string& usage,
                         const string& fileNotes,
                         bool visible,
                         const string& value,
                         const vector<string>& validValues)
  : Param(name, usage, fileNotes, visible),
    value_(value), original_(value), validValues_(validValues) {}
void StringParam::ThrowIfInvalid() const {
  if (!validValues_.empty() &&
      find(validValues_.begin(), validValues_.end(), value_) == validValues_.end()) {
    throw runtime_error("Invalid value for '" + name_ + "'; must be one of <" +
                        StringUtils::Join(validValues_, '|') + ">");
  }
}
string StringParam::GetType() const { return "string"; }
bool StringParam::IsDefault() const { return value_ == original_; }
bool StringParam::GetBool() const { return !value_.empty(); }
int StringParam::GetInt() const { return StringUtils::FromString<int>(value_); }
double StringParam::GetDouble() const { return StringUtils::FromString<double>(value_); }
string StringParam::GetString() const { return value_; }
void StringParam::Set(bool value) { value_ = value ? "true" : "false"; }
void StringParam::Set(int value) { value_ = StringUtils::ToString(value); }
void StringParam::Set(double value) { value_ = StringUtils::ToString(value); }
void StringParam::Set(const string& value) { value_ = value; }
//
// ArgParam
//
ArgParam::ArgParam(const string& name, const string& usage)
  : Param(name, usage, "", false), values_(vector<string>()) {}
string ArgParam::GetType() const { return "argument"; }
bool ArgParam::IsArgument() const { return true; }
bool ArgParam::IsDefault() const { return false; }
bool ArgParam::GetBool() const { return BoolParam::FromString(GetString()); }
int ArgParam::GetInt() const { return StringUtils::FromString<int>(GetString()); }
double ArgParam::GetDouble() const { return StringUtils::FromString<double>(GetString()); }
string ArgParam::GetString() const {
  if (values_.empty()) {
    throw runtime_error("No value for argument '" + name_ + "'");
  }
  return values_.front();
}
const vector<string>& ArgParam::GetStrings() const { return values_; }
void ArgParam::AddValue(const string& value) { values_.push_back(value); }

