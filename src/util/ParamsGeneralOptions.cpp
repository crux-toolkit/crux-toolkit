#include "Params.h"
#include "AminoAcidUtil.h"
#include "model/Peptide.h"
#include "model/objects.h"
#include "parameter.h"
#include "StringUtils.h"

#include <algorithm>

using namespace std;

void Params::InitGeneralOptions() {
  /* *** Initialize Options (command line and param file) *** */

  /* options for all executables */
  InitIntParam("verbosity", 30, 0, 100,
    "Specify the verbosity of the current processes. Each level prints the following "
    "messages, including all those at lower verbosity levels: 0-fatal errors, 10-non-"
    "fatal errors, 20-warnings, 30-information on the progress of execution, 40-more "
    "progress information, 50-debug info, 60-detailed debug info.",
    "Available for all crux programs.", true);
  InitStringParam("parameter-file", "",
    "A file containing parameters. [[html: See the "
    "<a href=\"../file-formats/parameter-file.html\">parameter documentation</a> page for details.]]",
    "Available for all crux programs. Any options specified on the "
    "command line will override values in the parameter file.", true);
  InitBoolParam("overwrite", false,
    "Replace existing files if true or fail when trying to overwrite a file if false.",
    "Available for all crux programs.  Applies to parameter file "
    "as well as index, search, and analysis output files.", true);
  /* generate_peptide parameters  */
  InitIntParam("min-length", 6, 1, MAX_PEPTIDE_LENGTH,
    "The minimum length of peptides to consider.",
    "Used from the command line or parameter file by "
    "crux-generate-peptides, and crux tide-index.", true);
  InitIntParam("max-length", 50, 1, MAX_PEPTIDE_LENGTH,
    "The maximum length of peptides to consider.",
    "Available from command line or parameter file for "
    "crux-generate-peptides, crux tide-index. ", true);
  InitDoubleParam("min-mass", 200, 0, BILLION,
    "The minimum mass (in Da) of peptides to consider.",
    "Available from command line or parameter file for "
    "crux-generate-peptides and crux tide-index. ", true);
  InitDoubleParam("max-mass", 7200, 1, BILLION,
    "The maximum mass (in Da) of peptides to consider.",
    "Available from command line or parameter file for "
    "crux-generate-peptides and crux tide-index. ", true);
  InitIntParam("min-peaks", 20, 0, BILLION,
    "The minimum number of peaks a spectrum must have for it to be searched.",
    "Available for tide-search.", true);
  InitStringParam("enzyme", "trypsin", "no-enzyme|trypsin|trypsin/p|chymotrypsin|"
    "elastase|clostripain|cyanogen-bromide|iodosobenzoate|proline-endopeptidase|"
    "staph-protease|asp-n|lys-c|lys-n|arg-c|glu-c|pepsin-a|"
    "elastase-trypsin-chymotrypsin|lysarginase|custom-enzyme",
    "Specify the enzyme used to digest the proteins in silico. Available enzymes "
    "(with the corresponding digestion rules indicated in parentheses) include "
    "no-enzyme ([X]|[X]), trypsin ([RK]|{P}), trypsin/p ([RK]|[]), chymotrypsin "
    "([FWYL]|{P}), elastase ([ALIV]|{P}), clostripain ([R]|[]), cyanogen-bromide "
    "([M]|[]), iodosobenzoate ([W]|[]), proline-endopeptidase ([P]|[]), staph-protease "
    "([E]|[]), asp-n ([]|[D]), lys-c ([K]|{P}), lys-n ([]|[K]), arg-c ([R]|{P}), "
    "glu-c ([DE]|{P}), pepsin-a ([FL]|{P}), elastase-trypsin-chymotrypsin "
    "([ALIVKRWFY]|{P}), lysarginase ([]|[KR]). "
    "Specifying --enzyme no-enzyme yields a non-enzymatic digest. "
    "[[html:<strong>]]Warning:[[html:</strong>]] the resulting index may be quite large.",
    "Available for crux-generate-peptides and crux tide-index.", true);
  InitStringParam("custom-enzyme", "",
    "Specify rules for in silico digestion of protein sequences. Overrides the enzyme "
    "option. Two lists of residues are given enclosed in square brackets or curly "
    "braces and separated by a |. The first list contains residues required/prohibited "
    "before the cleavage site and the second list is residues after the cleavage site. "
    "If the residues are required for digestion, they are in square brackets, '[' and "
    "']'. If the residues prevent digestion, then they are enclosed in curly braces, "
    "'{' and '}'. Use X to indicate all residues. For example, trypsin cuts after R or "
    "K but not before P which is represented as [RK]|{P}. AspN cuts after any residue "
    "but only before D which is represented as [X]|[D]. "
    "To prevent the sequences from being digested at all, use {X}|{X}.",
    "", true);
  InitDoubleParam("deisotope", 0, 0, 1000,
    "Perform a simple deisotoping operation across each MS2 spectrum. For each peak in an MS2 spectrum, consider lower m/z peaks. "
    "If the current peak occurs where an expected peak would lie for any charge state "
    "less than the charge state of the precursor, within mass tolerance, and if the "
    "current peak is of lower abundance, then the peak is removed.  The value of this "
    "parameter is the mass tolerance, in units of parts-per-million.  If set to 0, no "
    "deisotoping is performed.",
    "Available for tide-search.", true);
  InitStringParam("digestion", "full-digest",
    "full-digest|partial-digest|non-specific-digest",
    "Specify whether every peptide in the database must have two enzymatic termini "
    "(full-digest) or if peptides with only one enzymatic terminus are also included "
    "(partial-digest).",
    "Available for crux-generate-peptides and crux tide-index.", true);
  InitIntParam("missed-cleavages", 0, 0, 500,
    "Maximum number of missed cleavages per peptide to allow in enzymatic digestion.",
    "Available from command line or parameter file for crux-generate-peptides. "
    "When used with enzyme=<trypsin|elastase|chymotrypsin> "
    "includes peptides containing one or more potential cleavage sites.", true);
  InitDoubleParam("precursor-window", 50.0, 0, BILLION,
    "Tolerance used for matching peptides to spectra. Peptides must be within +/- "
    "'precursor-window' of the spectrum value. The precursor window units depend upon "
    "precursor-window-type.",
    "Available for tide-search and crux-generate-peptides.", true);
  InitStringParam("precursor-window-type", "ppm", "mass|mz|ppm",
    "Specify the units for the window that is used to select peptides around the precursor "
    "mass location (mass, mz, ppm). The magnitude of the window is defined by the precursor-"
    "window option, and candidate peptides must fall within this window. For the mass window-"
    "type, the spectrum precursor m+h value is converted to mass, and the window is defined "
    "as that mass +/- precursor-window. If the m+h value is not available, then the mass is "
    "calculated from the precursor m/z and provided charge. The peptide mass is computed as "
    "the sum of the monoisotopic amino acid masses plus 18 Da for the terminal OH group. The mz "
    "window-type calculates the window as spectrum precursor m/z +/- precursor-window and "
    "then converts the resulting m/z range to the peptide mass range using the precursor "
    "charge. For the parts-per-million (ppm) window-type, the spectrum mass is calculated as "
    "in the mass type. The lower bound of the mass window is then defined as the spectrum "
    "mass * (1.0 + (precursor-window / 1000000)) and the upper bound is defined as spectrum "
    "mass * (1.0 - (precursor-window / 1000000)).",
    "Available for tide-search.", true);
  InitStringParam("auto-precursor-window", "false", "false|warn|fail",
    "Automatically estimate optimal value for the precursor-window parameter "
    "from the spectra themselves. false=no estimation, warn=try to estimate "
    "but use the default value in case of failure, fail=try to estimate and "
    "quit in case of failure.",
    "Available for tide-search.", true);
  InitStringParam("spectrum-parser", "pwiz", "pwiz|mstoolkit",
    "Specify the parser to use for reading in MS/MS spectra.[[html: The default, "
    "ProteoWizard parser can read the MS/MS file formats listed <a href=\""
    "http://proteowizard.sourceforge.net/formats.shtml\">here</a>. The alternative is "
    "<a href=\"../mstoolkit.html\">MSToolkit parser</a>. "
    "If the ProteoWizard parser fails to read your files properly, you may want to try the "
    "MSToolkit parser instead.]]",
    "Available for tide-search.", true);
  InitBoolParam("use-z-line", true,
    "Specify whether, when parsing an MS2 spectrum file, Crux obtains the "
    "precursor mass information from the \"S\" line or the \"Z\" line. ",
    "Available when spectrum-parser = pwiz.", true);
  InitStringParam("keep-terminal-aminos", "NC", "N|C|NC|none",
    "When creating decoy peptides using decoy-format=shuffle or decoy-format="
    "peptide-reverse, this option specifies whether the N-terminal and "
    "C-terminal amino acids are kept in place or allowed to be shuffled or "
    "reversed. For a target peptide \"EAMPK\" with decoy-format=peptide-reverse, setting "
    "keep-terminal-aminos to \"NC\" will yield \"EPMAK\"; setting it to \"C\" will yield "
    "\"PMAEK\"; setting it to \"N\" will yield \"EKPMA\"; and setting it to \"none\" will "
    "yield \"KPMAE\".",
    "Available for tide-index.", true);
  InitBoolParam("peptide-list", false,
    "Create in the output directory a text file listing of all the peptides in the "
    "database, along with their corresponding decoy peptides, neutral masses and proteins, one per line.",
    "Available for tide-index.", true);
  // print-processed-spectra option
  InitStringParam("stop-after", "xcorr", "remove-precursor|square-root|"
    "remove-grass|ten-bin|xcorr",
    "Stop after the specified pre-processing step.",
    "Available for print-processed-spectra.", true);
  InitStringParam("output-units", "bin", "mz|bin",
    "Specify the output units for processed spectra.",
    "Available for print-processed-spectra", true);
  /* more generate_peptide parameters */
  InitBoolParam("sqt-output", false,
    "Outputs an SQT results file to the output directory. Note that if sqt-output is "
    "enabled, then compute-sp is automatically enabled and cannot be overridden.",
    "Available for tide-search.", true);
  InitBoolParam("mzid-output", false,
    "Output an mzIdentML results file to the output directory.",
    "Available for tide-search.", true);
  InitBoolParam("pin-output", false,
    "Output a Percolator input (PIN) file to the output directory.",
    "Available for tide-search.", true);
  InitBoolParam("mztab-output", false,
    "Output results in mzTab file to the output directory.",
    "Available for tide-search.", true);    
  InitBoolParam("pout-output", false,
    "Output a Percolator [[html:<a href=\""
    "https://github.com/percolator/percolator/blob/master/src/xml/percolator_out.xsd\">]]"
    "pout.xml[[html:</a>]] format results file to the output directory.",
    "Available for percolator.", true);
  InitBoolParam("pepxml-output", false,
    "Output a pepXML results file to the output directory.",
    "Available for tide-search, percolator.", true);
  InitBoolParam("txt-output", true,
    "Output a tab-delimited results file to the output directory.",
    "Available for tide-search, percolator.", true);
  InitStringParam("prelim-score-type", "sp", "sp|xcorr",
    "Initial scoring (sp, xcorr).",
    "The score applied to all possible psms for a given spectrum. Typically "
    "used to filter out the most plausible for further scoring.", false);
  InitStringParam("score-type", "xcorr", "xcorr|sp|xcorr-pvalue|sp-pvalue",
    "The primary scoring method to use (xcorr, sp, xcorr-pvalue, sp-pvalue).",
    "Primary scoring is typically done on a subset (see max-rank-preliminary) of all "
    "possible psms for each spectrum. Default is the SEQUEST-style xcorr. "
    "Crux also offers a p-value calculation for each psm based on xcorr "
    "or sp (xcorr-pvalue, sp-pvalue).", false);
  InitBoolParam("compute-sp", false,
    "Compute the preliminary score Sp for all candidate peptides. Report this score in the "
    "output, along with the corresponding rank, the number of matched ions and the total "
    "number of ions. This option is recommended if results are to be analyzed by Percolator."
    "If sqt-output is enabled, then compute-sp is automatically enabled and "
    "cannot be overridden. Note that the Sp computation requires re-processing each "
    "observed spectrum, so turning on this switch involves significant computational overhead.",
    "Available for tide-search.", true);
  InitStringParam("scan-number", "",
    "A single scan number or a range of numbers to be searched. Range should be "
    "specified as 'first-last' which will include scans 'first' and 'last'.",
    "The search range x-y is inclusive of x and y.", true);
  /* N.B. Use NaN to indicate that no user preference was specified.
   * In this case, the default value depends on the mass type.
   * S.M. Also prevent a width of 0.                                */
  InitDoubleParam("mz-bin-width", 0.02, 1e-4, BILLION,
    "Before calculation of the XCorr score, the m/z axes of the observed and theoretical "
    "spectra are discretized. This parameter specifies the size of each bin. The exact "
    "formula for computing the discretized m/z value is floor((x/mz-bin-width) + 1.0 - mz-bin-offset), where x is the observed m/z "
    "value. For low resolution ion trap ms/ms data 1.0005079 and for high resolution ms/ms "
    "0.02 is recommended.",
    "Available for tide-search.", true);
  InitDoubleParam("mz-bin-offset", 0.40, 0.0, 1.0,
    "In the discretization of the m/z axes of the observed and theoretical spectra, this "
    "parameter specifies the location of the left edge of the first bin, relative to "
    "mass = 0 (i.e., mz-bin-offset = 0.xx means the left edge of the first bin will be "
    "located at +0.xx Da).",
    "Available for tide-search.", true);
  InitStringParam("auto-mz-bin-width", "false", "false|warn|fail",
    "Automatically estimate optimal value for the mz-bin-width parameter "
    "from the spectra themselves. false=no estimation, warn=try to estimate "
    "but use the default value in case of failure, fail=try to estimate and "
    "quit in case of failure.",
    "Available for tide-search.", true);
  InitBoolParam("auto-modifications", false,
    "Automatically infer modifications from the spectra themselves.",
    "Available for tide-index.", true);
  InitStringParam("auto-modifications-spectra", "",
    "Specify the spectra file to be used for modification inference when the "
    "auto-modifications option is enabled. Multiple files may be separated by commas.",
    "Available for tide-index.", true);
  InitBoolParam("use-flanking-peaks", false,
    "Include flanking peaks around singly charged b and y theoretical ions. Each flanking "
    "peak occurs in the adjacent m/z bin and has half the intensity of the primary peak.",
    "Available for the tide-search command.", true);
  InitDoubleParam("spectrum-min-mz", 0.0, 0, BILLION,
    "The lowest spectrum m/z to search in the ms2 file.",
    "Available for tide-search.", true);
  InitDoubleParam("spectrum-max-mz", BILLION, 1, BILLION,
    "The highest spectrum m/z to search in the ms2 file.",
    "Available for tide-search.", true);
  InitStringParam("spectrum-charge", "all", "1|2|3|all",
    "The spectrum charges to search. With 'all' every spectrum will be searched and "
    "spectra with multiple charge states will be searched once at each charge state. "
    "With 1, 2, or 3 only spectra with that charge state will be searched.",
    "Used by tide-search.", true);
  InitStringParam("fileroot", "",
    "The fileroot string will be added as a prefix to all output file names.",
    "Available for all commands that produce an output directory.", true);
  InitStringParam("output-dir", "crux-output",
    "The name of the directory where output files will be created.",
    "Available for most commands.", true);
  InitStringParam("spectrum-outdir", "--output-dir value", "The name of the directory where the converted "
                  "spectrum files will be stored.", "Available for spectrum-converter", true);
  InitStringParam("temp-dir", "",
    "The name of the directory where temporary files will be created. If this "
    "parameter is blank, then the system temporary directory will be used",
    "Available for tide-index.", true);
  InitIntParam("memory-limit", 4, 1, BILLION, 
    "The maximum amount of memory (i.e., RAM), in GB, to be used by tide-index.",
    "Available for tide-index.", true);
  // coder options regarding decoys
  InitIntParam("num-decoy-files", 1, 0, 10,
    "Replaces number-decoy-set.  Determined by decoy-location"
    " and num-decoys-per-target",
    "", false);
  InitIntParam("num-decoys-per-target", 1, 1, BILLION,
    "The number of decoys to generate per target. When set to a value n, then "
    "with concat=F tide-search will output one target and n decoys. The "
    "resulting files can be used to run the \"average target-decoy "
    "competition\" method in assign-confidence. This parameter only applies "
    "when decoy-format=shuffle and should always be used in combination with "
    "allow-dups=T.",
    "Available for tide-index.", true);
  InitBoolParam("decoy-p-values", false,
    "Store all decoy p-values in a file",
    "", false);
  InitIntParam("top-match", 5, 1, BILLION,
    "Specify the number of matches to report for each spectrum.",
    "Available for tide-search and percolator", true);
  InitIntParam("top-match-in", 0, 0, BILLION,
    "Specify the maximum rank to allow when parsing results files. Matches with "
    "ranks higher than this value will be ignored (a value of zero allows matches with any rank).",
    "", true);
  InitStringParam("seed", "1",
    "When given a unsigned integer value seeds the random number generator with that value. "
    "When given the string \"time\" seeds the random number generator with the system time.",
    "Available for all Crux commands.", true);
  InitStringParam("fragment-mass", "mono", "average|mono",
    "Specify which isotopes to use in calculating fragment ion mass.",
    "Used by crux-predict-peptide-ions.", true);
  InitStringParam("isotopic-mass", "mono", "average|mono",
    "Specify the type of isotopic masses to use when calculating the peptide mass.",
    "Used from command line or parameter file by crux-generate-peptides.", true);
  InitIntParam("min-mods", 0, 0, MAX_PEPTIDE_LENGTH,
    "The minimum number of modifications that can be applied to a single "
    "peptide.",
    "Available for tide-index.", true);
  InitIntParam("max-mods", MAX_PEPTIDE_LENGTH, 0, MAX_PEPTIDE_LENGTH,
    "The maximum number of modifications that can be applied to a single "
    "peptide.",
    "Available for tide-index.", true);
  InitIntParam("max-aas-modified", MAX_PEPTIDE_LENGTH, 0, MAX_PEPTIDE_LENGTH,
    "The maximum number of modified amino acids that can appear in one "
    "peptide.  Each aa can be modified multiple times.",
    "", true);
  InitStringParam("mod-mass-format", "mod-only", "mod-only|total|separate",
    "Specify how sequence modifications are reported in various output files. Each "
    "modification is reported as a number enclosed in square braces following the "
    "modified residue; however, the number may correspond to one of three different "
    "masses: (1) 'mod-only' reports the value of the mass shift induced by the "
    "modification; (2) 'total' reports the mass of the residue with the modification "
    "(residue mass plus modification mass); (3) 'separate' is the same as 'mod-only', "
    "but multiple modifications to a single amino acid are reported as a "
    "comma-separated list of values. For example, suppose amino acid D has an "
    "unmodified mass of 115 as well as two moifications of masses +14 and +2. In this "
    "case, the amino acid would be reported as D[16] with 'mod-only', D[131] with 'total', "
    "and D[14,2] with 'separate'.",
    "Available for generate-peptides.", true);
  InitIntParam("mod-precision", 4, 0, 20,//arbitrary
    "Set the precision for modifications as written to .txt files.",
    "Also changes mods written to parameter file. By default, this "
    "value is set equal to the maximum modification precision in the "
    "specification of modifications.  Available for "
    "tide-index, tide-search and generate-peptides.",
    true);

  InitIntParam("precision", 8, 1, 100, //max is arbitrary
    "Set the precision for scores written to sqt and text files.",
    "Available for all commands.", true);
  InitIntParam("mass-precision", 4, 1, 100, // max is arbitrary
    "Set the precision for masses and m/z written to sqt and text files.",
    "Available for all commands.", true);
  InitIntParam("print-search-progress", 1000, 0, BILLION,
    "Show search progress by printing every n spectra searched. Set to 0 to show no "
    "search progress.",
    "Available for tide-search", true);
  // Sp scoring params
  InitDoubleParam("max-mz", 4000, 0, BILLION,
    "Used in scoring sp.",
    "Hide from users", false);
  InitDoubleParam("fraction-top-scores-to-fit", 0.55, 0, 1,
    "The fraction of psms per spectrum to use for estimating the "
    "score distribution for calculating p-values. "
    "Not compatible with 'number-top-scores-to-fig'.",
    "For developers/research only.", false);
  /* analyze-matches options */
  InitStringParam("algorithm", "percolator", "percolator|curve-fit|none",
    "The analysis algorithm to use (percolator, curve-fit, none).",
    "Available only for crux-analyze-matches.  Using 'percolator' will "
    "assign a q-value to the top-ranking psm for each spectrum based on "
    "the decoy searches.  Using 'curve-fit' will assign a q-value to same "
    "using the p-values calculated with score-type=<xcorr-pvalue|"
    "sq-pvalue>.  Incorrect combinations of score-type and algorithm cause"
    " undefined behavior. Using 'none' will turn the binary .csm files "
    "into text.", false);
  // **** percolator options. ****
  InitStringParam("search-input", "auto", "auto|separate|concatenated",
    "Specify the type of target-decoy search. Using 'auto', percolator attempts "
    "to detect the search type automatically.  Using 'separate' specifies two searches: "
    "one against target and one against decoy protein db. Using 'concatenated' "
    "specifies a single search on concatenated target-decoy protein db.",
    "Available for percolator", true);
  InitStringParam("percolator-seed", "1",
    "When given a unsigned integer value seeds the random number generator with that value. "
    "When given the string \"time\" seeds the random number generator with the system time.",
    "Available for all percolator", true);
InitStringParam("protein-name-separator", ",",
    "Determines the character to separate the protein IDs in the tab-delimited output format ",
    "Available for all percolator", true);
  InitBoolParam("feature-file-out", false,
    "Output the computed features in [[html:<a href=\"../file-formats/features.html\">]]"
    "tab-delimited Percolator input (.pin) format[[html:</a>]]. The features will be "
    "normalized, using either unit norm or standard deviation normalization (depending "
    "upon the value of the unit-norm option).",
    "Available for percolator.", true);
  InitBoolParam("decoy-xml-output", false,
    "Include decoys (PSMs, peptides, and/or proteins) in the XML output.",
    "Available for percolator", true);
  InitStringParam("decoy-prefix", "decoy_",
    "Specifies the prefix of the protein names that indicate a decoy.",
    "Available for tide-index and percolator", true);
  InitBoolParam("no-terminate", false,
    "Do not stop execution when encountering questionable SVM inputs or results. \"percolator.weights.txt\".",
    "Available for percolator", true); 
  InitBoolParam("output-weights", false,
    "Output final weights to a file named \"percolator.weights.txt\".",
    "Available for percolator", true);
  InitStringParam("init-weights", "",
    "Read the unnormalized initial weights from the third line of the given "
    "file. This can be the output of the --output-weights option from a "
    "previous Percolator analysis. Note that the weights must be in the same "
    "order as features in the PSM input file(s)",
    "Available for percolator", true);
  InitBoolParam("static", false,
    "Use the provided initial weights as a static model. If used, the "
    "--init-weights option must be specified.",
    "Available for percolator", true);
  InitIntParam("subset-max-train", 0,
    "Only train Percolator on a subset of PSMs, and use the resulting score "
    "vector to evaluate the other PSMs. Recommended when analyzing huge numbers "
    "(>1 million) of PSMs. When set to 0, all PSMs are used for training as "
    "normal.",
    "Available for percolator", true);
  InitDoubleParam("c-pos", 0.00,
    "Penalty for mistakes made on positive examples. If this value is set to 0, "
    "then it is set via cross validation over the values {0.1, 1, 10}, selecting the "
    "value that yields the largest number of PSMs identified at the q-value threshold "
    "set via the --test-fdr parameter.",
    "Available for percolator", true);
  InitDoubleParam("c-neg", 0.0, 0.0, 0.90,
    "Penalty for mistake made on negative examples. If not specified, then "
    "this value is set by cross validation over {0.1, 1, 10}.",
    "Available for percolator", true);
  InitDoubleParam("train-fdr", 0.01, 0, BILLION,
    "False discovery rate threshold to define positive examples in training.",
    "Available for percolator", true);
  InitDoubleParam("test-fdr", 0.01, 0.0, 1.0,
    "False discovery rate threshold used in selecting hyperparameters during internal "
    "cross-validation and for reporting the final results.",
    "Available for percolator.", true);
  InitDoubleParam("fido-fast-gridsearch", 0.0, 0.0, 1.0,
    "Apply the specified threshold to PSM, peptide and protein probabilities to "
    "obtain a faster estimate of the alpha, beta and gamma parameters.",
    "Available for percolator.", true);
  InitBoolParam("fido-no-split-large-components", false,
    "Do not approximate the posterior distribution by allowing large graph "
    "components to be split into subgraphs. The splitting is done by "
    "duplicating peptides with low probabilities. Splitting continues "
    "until the number of possible configurations of each subgraph is "
    "below 2^18",
    "Available for percolator", true);
  InitDoubleParam("fido-protein-truncation-threshold", 0.01, 0.0, 1.0,
    "To speed up inference, proteins for which none of the associated "
    "peptides has a probability exceeding the specified threshold will "
    "be assigned probability = 0.",
    "Available for percolator", true);
  InitBoolParam("tdc", true,
    "Use target-decoy competition to assign q-values and PEPs. When set to F, "
    "the mix-max method, which estimates the proportion pi0 of incorrect target "
    "PSMs, is used instead.",
    "Available for percolator", true);
  InitIntParam("maxiter", 10, 0, 100000000,
    "Maximum number of iterations for training.",
    "Available for percolator", true);
  InitBoolParam("quick-validation", false,
    "Quicker execution by reduced internal cross-validation.",
    "Available for percolator", true);
  InitStringParam("default-direction", "",
    "In its initial round of training, Percolator uses one feature to induce a ranking "
    "of PSMs. By default, Percolator will select the feature that produces the largest "
    "set of target PSMs at a specified FDR threshold (cf. --train-fdr). This option "
    "allows the user to specify which feature is used for the initial ranking, using the "
    "name as a string[[html: from <a href=\"../file-formats/features.html\">this table</a>]]. The name "
    "can be preceded by a hyphen (e.g. \"-XCorr\") to indicate that a lower value is "
    "better.",
    "Available for percolator", true);
  InitBoolParam("unitnorm", false,
    "Use unit normalization (i.e., linearly rescale each PSM's feature vector to have a "
    "Euclidean length of 1), instead of standard deviation normalization.",
    "Available for percolator.", true);
  InitBoolParam("test-each-iteration", false,
    "Measure performance on test set each iteration.",
    "Available for percolator.", true);
  InitStringParam("picked-protein", "",
    "Use the picked protein-level FDR to infer protein probabilities, provide the "
    "fasta file as the argument to this flag.",
    "Available for percolator", true);
  InitStringParam("protein-enzyme", "trypsin", "no_enzyme|elastase|pepsin|proteinasek|"
    "thermolysin|trypsinp|chymotrypsin|lys-n|lys-c|arg-c|asp-n|glu-c|lysarginase|trypsin",
    "Type of enzyme",
    "Available for percolator", true);
  InitBoolParam("protein-report-fragments", false,
    "By default, if the peptides associated with protein A are a proper subset "
    "of the peptides associated with protein B, then protein A is eliminated and "
    "all the peptides are considered as evidence for protein B. Note that this "
    "filtering is done based on the complete set of peptides in the database, not "
    "based on the identified peptides in the search results. Alternatively, if this "
    "option is set and if all of the identified peptides associated with protein B "
    "are also associated with protein A, then Percolator will report a comma-"
    "separated list of protein IDs, where the full-length protein B is first in the "
    "list and the fragment protein A is listed second. Not available for Fido.",
    "Available for percolator", true);
  InitBoolParam("protein-report-duplicates", false,
    "If multiple database proteins contain exactly the same set of peptides, then "
    "Percolator will randomly discard all but one of the proteins. If this option "
    "is set, then the IDs of these duplicated proteins will be reported as a comma-"
    "separated list. Not available for Fido.",
    "Available for percolator", true);
  InitBoolParam("protein", false,
    "Use the Fido algorithm to infer protein probabilities. Must be true to use any of the Fido options.",
    "Available for percolator", true);
  InitDoubleParam("fido-alpha", 0.0, 0.0, 1.0,
    "Specify the probability with which a present protein emits an associated peptide. "
    "Set by grid search (see --fido-gridsearch-depth parameter) if not specified.",
    "Available for percolator if --protein T is set.", true);
  InitDoubleParam("fido-beta", 0.0, 0.0, 10.0,
    "Specify the probability of the creation of a peptide from noise. Set by grid "
    "search (see --fido-gridsearch-depth parameter) if not specified.",
    "Available for percolator if --protein T is set.", true);
  InitDoubleParam("fido-gamma", 0.0, 0.0, 10.0,
    "Specify the prior probability that a protein is present in the sample. Set by grid "
    "search (see --fido-gridsearch-depth parameter) if not specified.",
    "Available for percolator if --protein T is set.", true);
  InitBoolParam("fido-empirical-protein-q", false,
    "Estimate empirical p-values and q-values for proteins using target-decoy analysis.",
    "Available for percolator if --protein T is set.", true);
  InitIntParam("fido-gridsearch-depth", 0, 0, 2,
    "Set depth of the grid search for alpha, beta and gamma estimation.[[html: The values "
    "considered, for each possible value of the --fido-gridsearch-depth parameter, are as follows:<ul>"
    "<li>0: alpha = {0.01, 0.04, 0.09, 0.16, 0.25, 0.36, 0.5}; beta = {0.0, 0.01, 0.15, "
    "0.025, 0.035, 0.05, 0.1}; gamma = {0.1, 0.25, 0.5, 0.75}.</li><li>1: alpha = {0.01, "
    "0.04, 0.09, 0.16, 0.25, 0.36}; beta = {0.0, 0.01, 0.15, 0.025, 0.035, 0.05}; gamma = "
    "{0.1, 0.25, 0.5}.</li><li>2: alpha = {0.01, 0.04, 0.16, 0.25, 0.36}; beta = {0.0, "
    "0.01, 0.15, 0.030, 0.05}; gamma = {0.1, 0.5}.</li><li>3: alpha = {0.01, 0.04, 0.16, "
    "0.25, 0.36}; beta = {0.0, 0.01, 0.15, 0.030, 0.05}; gamma = {0.5}.</li></ul>]]",
    "Available for percolator if --protein T is set.", true);
  InitDoubleParam("fido-gridsearch-mse-threshold", 0.05, 0, 1,
    "Q-value threshold that will be used in the computation of the MSE and ROC AUC "
    "score in the grid search.",
    "Available for percolator if --protein T is set.", true);
  InitBoolParam("override", false,
    "By default, Percolator will examine the learned weights for each feature, and if "
    "the weight appears to be problematic, then percolator will discard the learned "
    "weights and instead employ a previously trained, static score vector. This switch "
    "allows this error checking to be overriden.",
    "Available for percolator.", true);
  InitBoolParam("klammer", false,
    "Use retention time features calculated as in \"Improving tandem mass spectrum "
    "identification using peptide retention time prediction across diverse chromatography "
    "conditions\" by Klammer AA, Yi X, MacCoss MJ and Noble WS. ([[html:<em>]]Analytical "
    "Chemistry[[html:</em>]]. 2007 Aug 15;79(16):6111-8.).",
    "Available for percolator", true);
  InitBoolParam("only-psms", false,
    "Report results only at the PSM level.  This flag causes Percolator to skip the "
    "step that selects the top-scoring PSM per peptide; hence, peptide-level results "
    "are left out and only PSM-level results are reported.",
    "Available for percolator", true);
  InitBoolParam("train-best-positive", false,
    "Enforce that, for each spectrum, at most one PSM is included in the "
    "positive set during each training iteration. Note that if the user only "
    "provides one PSM per spectrum, then this option will have no effect.",
    "Available for percolator", true);
  InitDoubleParam("spectral-counting-fdr", 0, 0, 1,
    "Report the number of unique PSMs and total (including shared peptides) "
    "PSMs as two extra columns in the protein tab-delimited output.",
    "Available for percolator", true);
  // **** Tide arguments ****
  InitArgParam("spectrum records file",
    "A spectrum records file generated by a previous run of crux tide-search "
    "using the store-spectra parameter.");
  InitArgParam("tide spectra file",
    "The name of one or more files from which to parse the fragmentation spectra, in any "
    "of the file formats supported by ProteoWizard. Alternatively, the argument "
    "may be one or more binary spectrum files produced by a previous run of crux "
    "tide-search using the store-spectra parameter. Multiple files can be included "
    "on the command line (space delimited), prior to the name of the database.");
  InitArgParam("tide database",
    "Either a FASTA file or a directory containing a database index created by a previous "
    "run of crux tide-index.");
  // **** Tide options ****
  InitStringParam("decoy-format", "shuffle", "none|shuffle|peptide-reverse",
    "Include a decoy version of every peptide by shuffling or reversing the "
    "target sequence or protein. In shuffle or peptide-reverse mode, each peptide is "
    "either reversed or shuffled, leaving the N-terminal and C-terminal amino acids in "
    "place. Note that peptides appear multiple times in the target database are only "
    "shuffled once. In peptide-reverse mode, palindromic peptides are shuffled. Also, if a "
    "shuffled peptide produces an overlap with the target or decoy database, then the "
    "peptide is re-shuffled up to 5 times. Note that, despite this repeated shuffling, "
    "homopolymers will appear in both the target and decoy database.",
    "Available for tide-index", true);
  InitStringParam("mods-spec", "C+57.02146",
    "[[nohtml:Expression for static and variable mass modifications to include. "
    "Specify a comma-separated list of modification sequences of the form: "
    "C+57.02146,2M+15.9949,1STY+79.966331,...]][[html:The general form of a modification "
    "specification has three components, as exemplified by <span style=\"color: red;\">1"
    "</span><span style=\"color: green;\">STY</span>+<span style=\"color: blue\">79.966331"
    "</span>.<br>The three components are: [<span style=\"color: red;\">max_per_peptide"
    "</span>]<span style=\"color: green;\">residues</span>[+/-]<span style-\"color: blue;\">"
    "mass_change</span><br>In the example, <span style=\"color: red;\">max_per_peptide"
    "</span> is <span style=\"color: red;\">1</span>, <span style=\"color: green;\">"
    "residues</span> are <span style=\"color: green;\">STY</span>, and "
    "<span style=\"color: blue;\">mass_change</span> is <span style=\"color: blue;\">"
    "+79.966331</span>. To specify a static modification, the number preceding the amino "
    "acid must be omitted; i.e., <span style=\"color: green;\">C</span>+<span "
    "style=\"color: blue;\">57.02146</span> specifies a static modification of 57.02146 "
    "Da to cysteine. Note that Tide allows at most one modification per amino "
    "acid.  Also, the default modification (C+57.02146) will be added to "
    "every mods-spec string unless an explicit C+0 is included.]]",
    "Available for tide-index.", true);
  InitStringParam("nterm-peptide-mods-spec", "",
    "[[nohtml:Specifies N-terminal static and variable mass modifications on peptides. "
    "Specify a comma-separated list of N-terminal modification sequences of the form: "
    "1E-18.0106,C-17.0265]][[html:Specify peptide n-terminal modifications. Like "
    "--mods-spec, this specification has three components, but with a slightly different "
    "syntax. The <span style=\"color: red;\">max_per_peptide</span> can be either \"1\", "
    "in which case it defines a variable terminal modification, or missing, in which case "
    "the modification is static. The <span style=\"color: green;\">residues</span> field "
    "indicates which amino acids are subject to the modification, with the residue <span "
    "style=\"color: green;\">X</span> corresponding to any amino acid. Finally, <span "
    "style=\"color: blue;\">added_mass</span> is defined as before.]]",
    "Available for tide-index", true);
  InitStringParam("cterm-peptide-mods-spec", "",
    "[[nohtml:Specifies C-terminal static and variable mass modifications on peptides. "
    "Specify a comma-separated list of C-terminal modification sequences of the form: "
    "X+21.9819]][[html:Specify peptide c-terminal modifications. See "
    "nterm-peptide-mods-spec for syntax.]]",
    "Available for tide-index", true);
  InitStringParam("cterm-protein-mods-spec", "",
    "[[nohtml:Specifies static and variable mass modifications on proteins' C-terminal. ]]"
    "[[html:Specifies C-terminal static and variable mass modifications on proteins."
    "Mod specification syntax is the same as for peptide mods (see nterm-peptide-mods-spec option),"
    "but these mods are applied only to peptide C-terminals that are also protein terminals."
    "If variable modification are provided for both peptide and protein terminal, they will be "
    "applied one at a time. ",
    "Available for tide-index", true);
  InitStringParam("nterm-protein-mods-spec", "",
    "[[nohtml:Specifies static and variable mass modifications on proteins N-terminal.]] "
    "[[html:Same as cterm-protein-mods-spec, but for  the protein N-terminal.",
    "Available for tide-index", true);
  InitStringParam("store-spectra", "",
    "Specify the name of the file where the binarized fragmentation spectra "
    "will be stored. Subsequent runs of crux tide-search will execute more quickly if "
    "provided with the spectra in binary format. The filename is specified relative to "
    "the current working directory, not the Crux output directory (as specified by "
    "--output-dir). This option is not valid if multiple input spectrum files are given.",
    "Available for tide-search", true);
  InitBoolParam("exact-p-value", false,
    "Enable the calculation of exact p-values for the XCorr score[[html: as described in "
    "<a href=\"http://www.ncbi.nlm.nih.gov/pubmed/24895379\">this article</a>]]. Calculation "
    "of p-values increases the running time but increases the number of identifications at a "
    "fixed confidence threshold. The p-values will be reported in a new column with header "
    "\"exact p-value\", and the \"xcorr score\" column will be replaced with a \"refactored "
    "xcorr\" column. Note that, currently, p-values can only be computed when the "
    "mz-bin-width parameter is set to its default value. Variable and static mods are allowed "
    "on non-terminal residues in conjunction with p-value computation, but currently only "
    "static mods are allowed on the N-terminus, and no mods on the C-terminus.",
    "Available for tide-search", true);
  //Added for tailor score calibration method by AKF
  InitBoolParam("use-tailor-calibration", false,
    "Fast, but heuristic PSM score calibration[[html: as described in "
    "<a href=\"https://pubmed.ncbi.nlm.nih.gov/32175744/\">this article</a>]].",
    "Available for tide-search", true);    
  InitStringParam("store-index", "",
    "When providing a FASTA file as the index, the generated binary index will be stored at "
    "the given path. This option has no effect if a binary index is provided as the index.",
    "Available for tide-search", true);
  InitBoolParam("concat", false,
    "When set to T, target and decoy search results are reported in a single file, and only "
    "the top-scoring N matches (as specified via --top-match) are reported for each spectrum, "
    "irrespective of whether the matches involve target or decoy peptides.",
    "Available for tide-search", true);
  InitBoolParam("file-column", true,
    "Include the file column in tab-delimited output.",
    "Available for tide-search", true);
  // Same as remove_precursor_peak and remove_precursor tolerance in Comet
  InitBoolParam("remove-precursor-peak", false,
    "If true, all peaks around the precursor m/z will be removed, within a range "
    "specified by the --remove-precursor-tolerance option.",
    "Available for tide-search.", true);
  InitDoubleParam("remove-precursor-tolerance", 1.5, 0, BILLION,
    "This parameter specifies the tolerance (in Th) around each precursor m/z that is "
    "removed when the --remove-precursor-peak option is invoked.",
    "Available for print-processed spectra and tide-search.", true);
  InitBoolParam("clip-nterm-methionine", false,
    "When set to T, for each protein that begins with methionine, tide-index will "
    "put two copies of the leading peptide into the index, with and without the N-terminal "
    "methionine.",
    "Available for tide-index.", true);
  InitBoolParam("allow-dups", false,
    "Prevent duplicate peptides between the target and decoy databases. When set to \"F\", "
    "the program keeps all target and previously generated decoy peptides in memory. A shuffled "
    "decoy will be re-shuffled multiple times to avoid duplication. If a non-duplicated peptide "
    "cannot be generated, the decoy is skipped entirely. When set to \"T\", every decoy is added to "
    "the database without checking for duplication. This option reduces the memory requirements "
    "significantly.",
    "Available for tide-index.", true);
  InitBoolParam("use-neutral-loss-peaks", true,
    "Controls whether neutral loss ions are considered in the search. "
    "For XCorr, the loss of ammonia (NH3, 17.0086343 Da) is applied to singly "
    "charged b- and y-ions, and the loss of water (H2O; 18.0091422) is applied "
    "to b-ions. If the precursor charge is >=3, then a doubly-charged version of "
    "each ion is added. For XCorr p-value, three types of neutral losses are included. "
    "Loss of ammonia and water are applied to b- and y-ions, and a carbon monoxide "
    "loss (CO, 27.9949) is also applied to b-ions. Higher charge fragments are "
    "included for all possible charges less than the precursor charge. All "
    "neutral loss peaks have an intensity 1/10 of the primary peak. Neutral losses "
    "are not yet implemented for the res-ev score function.",
    "Available for tide-search.", true);
  InitIntParam("min-precursor-charge", 1, 1, BILLION,
    "The minimum charge state of a spectra to consider in search.",
    "Available for tide-search.", true);
  InitIntParam("max-precursor-charge", 5, 1, BILLION,
    "The maximum charge state of a spectra to consider in search.",
    "Available for tide-search.", true);
  InitBoolParam("peptide-centric-search", false,
    "Carries out a peptide-centric search. For each peptide the top-scoring spectra "
    "are reported, in contrast to the standard spectrum-centric search where the top-"
    "scoring peptides are reported. Note that in this case the \"xcorr rank\" column "
    "will contain the rank of the given spectrum with respect to the given candidate "
    "peptide, rather than vice versa (which is the default).",
    "Available for tide-search.", true);
  InitIntParam("elution-window-size", 0, 0, 10,
    "Size of the elution window used in smoothing score in DIA mode. "
    "Used only with peptide-centric-search if greater than 0. A score of a psms "
    "centred in the window is substituted by the geometric mean of the scores "
    "in the window. If windows size is even, then it is increased by 1.",
    "Available for tide-search.", false);
  InitBoolParam("skip-decoys", true,
    "Skips decoys when reading a Tide index.",
    "Available for read-tide-index", false);
  InitBoolParam("skip-preprocessing", false,
    "Skip preprocessing steps on spectra. Default = F.",
    "Available for tide-search", true);
  InitStringParam("score-function", "xcorr", "xcorr|combined-p-values",  //TDOD: implement later : |hyperscore|hyperscore-la
    "Function used for scoring PSMs. 'xcorr' is the original scoring function used by SEQUEST;"
    "`combined-p-values` combined (1) exact-p-value: a calibrated version of XCorr that uses "
    "dynamic programming and (2) residue-evidence-pvalue: a valibarated version of the  ResEV "
    "that considers pairs of peaks, rather than single peaks; ",
//    "`hyperscore` is the score function used in X!Tandem; `hyperscore-la` is a variant of the "  TODO: implement later
//    "hyperscore designed for open modification searching.",
    "Available for tide-search.", true);
  InitDoubleParam("fragment-tolerance", .02, 0, 2,
    "Mass tolerance (in Da) for scoring pairs of peaks when creating the residue evidence matrix. "
    "This parameter only makes sense when score-function is 'residue-evidence' or 'both'.",
    "Available for tide-search.", true);
  InitIntParam("evidence-granularity", 25, 1, 100,
    "This parameter controls the granularity of the entries in the dynamic programming matrix used in residue-evidence scoring."
    "Smaller values make the program run faster but give less accurate p-values; "
    "larger values make the program run more slowly but give more accurate p-values.",
    "Available for tide-search", true);
  InitStringParam("isotope-error", "",
                  "List of positive, non-zero integers.",
                  "Isotope errors to include. "
                  "Specify a comma-separated list of isotope errors of the form: "
                  "1,2,3,..."
                  "Available for tide-search", true);
  InitIntParam("num-threads", 1, 0, 64,
               "0=poll CPU to set num threads; else specify num threads directly.",
               "Available for tide-search tab-delimited files only.", true);
  InitBoolParam("brief-output", false,
    "Output in tab-delimited text only the file name, scan number, charge, score and peptide."
    "Incompatible with mzid-output=T, pin-output=T, pepxml-output=T or txt-output=F.",
    "Available for tide-search", true);
  InitBoolParam("override-charges", false,
    "If this is set to T, then all spectra are searched in all charge states from min-charge to max-charge. "
    "Otherwise, the default behavior is to search with all charge states only if a spectrum has no charge "
    "or charge=0.",
    "Available for tide-search", true);
}
