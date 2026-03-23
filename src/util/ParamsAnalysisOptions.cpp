#include "Params.h"
#include "AminoAcidUtil.h"
#include "model/Peptide.h"
#include "model/objects.h"
#include "parameter.h"
#include "StringUtils.h"

#include <algorithm>

using namespace std;

void Params::InitAnalysisOptions() {
  /* analyze-matches parameter options */
  InitArgParam("target input",
    "One or more files, each containing a collection of peptide-spectrum matches (PSMs) "
    "in [[html:<a href=\"../file-formats/txt-format.html\">]]tab-delimited text[[html:</a>]], [[html:<a "
    "href=\"http://tools.proteomecenter.org/wiki/index.php?title=Formats:pepXML\">]]PepXML"
    "[[html:</a>]], or [[html:<a href=\"http://www.psidev.info/mzidentml\">]]mzIdentML"
    "[[html:</a>]] format. In tab-delimited text format, only the specified score column "
    "is required. However if --estimation-method is tdc, then the columns \"scan\" and "
    "\"charge\" are required, as well as \"protein ID\" if the search was run with "
    "concat=F. Furthermore, if the --estimation-method is specified to peptide-level "
    "is set to T, then the column "
    "\"peptide\" must be included, and if --sidak is set to T, then the \"distinct "
    "matches/spectrum\" column must be included.[[html:<br>Note that multiple files can "
    "also be provided either on the command line or using the --list-of-files option.<br>"
    "Decoys can be provided in two ways: either as a separate file or embedded within the "
    "same file as the targets. Crux will first search the given file for decoys using a "
    "prefix (specified via --decoy-prefix) on the protein name. If no decoys are found, "
    "then Crux will search for decoys in a separate file. The decoy file name is constructed "
    "from the target file name by replacing \"target\" with \"decoy\". For example, if "
    "tide-search.target.txt is provided as input, then Crux will search for a corresponding "
    "file named \"tide-search.decoy.txt.\"<br>Note that if decoys are provided in a separate "
    "file, then assign-confidence will first carry out a target-decoy competition, "
    "identifying corresponding pairs of targets and decoys and eliminating the one with "
    "the worse score. In this case, the column/tag called \"delta_cn\" will be eliminated "
    "from the output.]]");
  InitDoubleParam("pi-zero", 1.0, 0, 1,
    "The estimated percent of target scores that are drawn from the "
    "null distribution.",
    "Used by assign-confidence and percolator", false);
  InitStringParam("estimation-method", "tdc", "mix-max|tdc|peptide-level",
    "Specify the method used to estimate q-values.  The mix-max procedure or target-decoy "
    "competition apply to PSMs. The peptide-level option eliminates any PSM for which there "
    "exists a better scoring PSM involving the same peptide, and then uses decoys to "
    "assign confidence estimates.",
    "Used by assign-confidence.", true);
  InitBoolParam("sidak", false,
    "Adjust the score using the Sidak adjustment and reports them in a new column in the "
    "output file. Note that this adjustment only makes sense if the given scores are "
    "p-values, and that it requires the presence of the \"distinct matches/spectrum\" "
    "feature for each PSM.",
    "Used by assign-confidence.", true);
  InitStringParam("score", "",
    "Specify the column (for tab-delimited input) or tag (for XML input) "
    "used as input to the q-value estimation procedure. If this parameter is unspecified, "
    "then the program searches for \"xcorr score\", \"evalue\" (comet), "
    "\"exact p-value\" score fields in this order in the input file. ",
    "Used by assign-confidence.", true);
  InitBoolParam("combine-charge-states", false,
    "Specify this parameter to T in order to combine charge states with peptide sequences"
    "in peptide-centric search. Works only if estimation-method = peptide-level.",
    "Used by assign-confidence.", true);
  InitBoolParam("combine-modified-peptides", false,
    "Specify this parameter to T in order to treat peptides carrying different or "
    "no modifications as being the same. Works only if estimation = peptide-level.",
    "Used by assign-confidence.", true);
  InitStringParam("percolator-intraset-features", "F",
    "Set a feature for percolator that in later versions is not an option.",
    "Shouldn't be variable; hide from user.", false);
  InitBoolParam("use-old-atdc", false,
                "Use the originally described version of aTDC, rather than the improved one.",
                "Used by assign-confidence.", false);
  /* Cascade-Search parameters */
  InitDoubleParam("q-value-threshold", 0.01, 0, 1.0,
    "The q-value threshold used by cascade search. Each spectrum identified in one search "
    "with q-value less than this threshold will be excluded from all subsequent searches. "
    "Note that the threshold is not applied to the final database in the cascade.",
    "Used by cascade-search.", true);
  InitArgParam("database-series",
    "A comma-separated list of databases, each generated by tide-index. "
    "Cascade-search will search the given spectra against these databases in the given order.");
  InitIntParam("cascade-termination", 20, 0, BILLION,
    "The minimum number of accepted PSMs required for cascade-search to continue to the "
    "next database in the given series",
    "Used by cascade-search.", false);
  /*Subtract-index parameters*/
  InitArgParam("tide index 1", "A peptide index produced using tide-index");
  InitArgParam("tide index 2", "A second peptide index, to be subtracted from the first index.");
  InitArgParam("output index", "A new peptide index containing all peptides that occur in the"
    "first index but not the second.");
//  InitArgParam("index name", "output tide index");
  // **** predict-peptide-ions options. ****
  InitStringParam("primary-ions", "by", "a|b|y|by|bya",
    "Predict the specified primary ion series. 'a' indicates a-ions only, 'b' indicates "
    "b-ions only, 'y' indicates y-ions only, 'by' indicates both b and y, 'bya' "
    "indicates b, y, and a.",
    "Only available for crux-predict-peptide-ions. Set automatically to "
    "'by' for searching.", true);
  InitBoolParam("precursor-ions", false,
    "Predict the precursor ions, and all associated ions (neutral losses, multiple "
    "charge states) consistent with the other specified options.",
    "Only available for crux-predict-peptide-ions.", true);
  InitIntParam("isotope", 0, 0, 2,
    "Predict the given number of isotope peaks (0|1|2).",
    "Only available for crux-predict-peptide-ion.  Automatically set to "
    "0 for Sp scoring and 1 for xcorr scoring.", true);
  InitBoolParam("flanking", false,
    "Predict flanking peaks for b- and y ions.",
    "Only available for crux-predict-peptide-ion.", true);
  InitStringParam("max-ion-charge", "peptide",
    "Predict theoretical ions up to max charge state (1, 2, ... ,6) or up to the charge state "
    "of the peptide (\"peptide\"). If the max-ion-charge is greater than the "
    "charge state of the peptide, then the maximum is the peptide charge. ",
    "Available for predict-peptide-ions. "
    "Set to 'peptide' for search.", true);
  InitIntParam("nh3", 0, -100, BILLION,
    "Include among the predicted peaks b/y ions with up to n losses of nh3. For example, "
    "for --nh3 2, predict a peak for each b- and y-ion with the loss of one nh3 group and "
    "predict a second peak for each b- and y-ion with the loss of two nh3 groups. These "
    "peaks will have 1 and 2, respectively, in the NH3 column of the output.",
    "Only available for crux-predict-peptide-ions.", true);
  InitIntParam("h2o", 0, -100, BILLION,
    "Include in the predicted peaks, b/y ions with the loss of 1 to n water molecules. See "
    "--nh3 for an example.",
    "Only available for crux-predict-peptide-ions.", true);
  // ***** spectral-counts aguments *****
  InitArgParam("input PSMs",
    "A PSM file in either tab delimited text format (as produced by percolator), or pepXML format.");
  // also uses "protein-database"
  // ***** spectral-counts options *****
  InitStringParam("protein-database", "",
    "The name of the file in FASTA format.",
    "Option for spectral-counts", true);
  InitStringParam("measure", "NSAF", "RAW|NSAF|dNSAF|SIN|EMPAI",
    "Type of analysis to make on the match results: "
    "(RAW|NSAF|dNSAF|SIN|EMPAI). With exception of the RAW metric, the database of "
    "sequences need to be provided using --protein-database.",
    "Available for spectral-counts.  RAW is raw counts, "
    "NSAF is Normalized Spectral Abundance Factor, "
    "dNSAF is Distributed Spectral Abundance Factor, "
    "SIN is Spectral Index Normalized and EMPAI is "
    "Exponentially Modified Protein Abundance Index", true);
  InitBoolParam("unique-mapping", false,
    "Ignore peptides that map to multiple proteins.",
    "Available for spectral-counts.", true);
  InitStringParam("quant-level", "protein", "protein|peptide",
    "Quantification at protein or peptide level.",
    "Available for spectral-counts and either NSAF and SIN.", true);
  InitStringParam("parsimony", "none", "none|simple|greedy",
    "Perform a parsimony analysis on the proteins, and report a "
    "\"parsimony rank\" column in the output file. This column contains "
    "integers indicating the protein's rank in a list sorted by spectral "
    "counts. If the parsimony analysis results in two proteins being merged, "
    "then their parsimony rank is the same. In such a case, the rank is "
    "assigned based on the largest spectral count of any protein in the merged "
    "meta-protein. The \"simple\" parsimony algorithm only merges two proteins "
    "A and B if the peptides identified in protein A are the same as or a "
    "subset of the peptides identified in protein B. The \"greedy\" parsimony "
    "algorithm does additional merging, by identifying the longest protein "
    "(i.e., the protein with the most peptides) that contains one or more "
    "shared peptides. The shared peptides are assigned to the identified "
    "protein and removed from any other proteins that contain them, and the "
    "process is then repeated. Note that, with this option, some proteins end "
    "up being assigned no peptides at all; these orphan proteins are not "
    "reported in the output.", 
    "Available for spectral-counts.", true);
  InitStringParam("threshold-type", "qvalue", "none|qvalue|custom",
    "Determines what type of threshold to use when filtering matches. none : read all "
    "matches, qvalue : use calculated q-value from percolator, custom : use "
    "--custom-threshold-name and --custom-threshold-min parameters.",
    "used for crux spectral-counts", true);
  InitDoubleParam("threshold", 0.01,
    "Only consider PSMs with a threshold value. By default, q-values "
    "are thresholded using a specified threshold value. This behavior can be "
    "changed using the --custom-threshold and --threshold-min "
    "parameters.",
    "Available for spectral-counts. All PSMs with higher (or lower) than "
    "this will be ignored.", true);
  InitStringParam("custom-threshold-name", "",
    "Specify which field to apply the threshold to. The direction of the threshold "
    "(<= or >=) is governed by --custom-threshold-min. By default, the threshold "
    "applies to the q-value, specified by \"percolator q-value\", "
    "\"decoy q-value (xcorr)\".",
    "Available for spectral-counts.", true);
  InitBoolParam("custom-threshold-min", true,
    "When selecting matches with a custom threshold, custom-threshold-min determines "
    "whether to filter matches with custom-threshold-name values that are greater-than or "
    "equal (F) or less-than or equal (T) than the threshold.",
    "Available for spectral-counts.", true);
  InitStringParam("input-ms2", "",
    "MS2 file corresponding to the psm file. Required to measure the SIN. Ignored for "
    "NSAF, dNSAF and EMPAI.",
    "Available for spectral-counts with measure=SIN.", true);
  InitBoolParam("mzid-use-pass-threshold", false,
    "Use mzid's passThreshold attribute to filter matches.",
    "Used when parsing mzIdentML files.", true);
  InitBoolParam("find-peptides", true,
                "Validate peptides by finding them in the given FASTA file.",
                "Only available for spectral-counts.", false);
  // ***** static mods *****
  for (char c = 'A'; c <= 'Z'; c++) {
    double deltaMass = (c != 'C') ? 0 : CYSTEINE_DEFAULT;
    bool visible = (c != 'B' && c != 'J' && c != 'O' && c != 'U' && c != 'X' && c != 'Z');
    InitDoubleParam(string(1, c), deltaMass,
      "Change the mass of all amino acids '" + string(1, c) + "' by the "
      "given amount.", "", visible);
  }
  /* psm-convert options */
  InitStringParam("input-format", "auto", "auto|tsv|sqt|pepxml|mzidentml",
    "Legal values are auto, tsv, sqt, pepxml or mzidentml format.",
    "option, for psm-convert", true);
  InitBoolParam("distinct-matches", true,
    "Whether matches/ion are distinct (as opposed to total).",
    "option, for psm-convert.", true);
  /* get-ms2-spectrum options */
  InitBoolParam("stats", false,
    "Rather than the spectrum, output summary statistics to standard output. Each statistic "
    "is placed on a separate line, in the format <name>:<value> (e.g. \"TIC:1000.0\")."
    "[[html:<br>The following statistics are reported for the entire spectrum:<ul><li>"
    "Precursor m/z</li><li>Total Ion Current</li><li>Base Peak Intensity</li><li>Number of "
    "peaks</li><li>Minimum m/z</li><li>Maximum m/z</li></ul>In addition, for each possible "
    "spectrum charge state, the following statistics are reported:<ul><li>Charge state</li>"
    "<li>Neutral mass</li><li>Charged mass</li><li>M+H+ mass</li></ul>]]",
    "Available only for crux-get-ms2-spectrum.  Does not affect contents "
    "of the output file.", true);

  /* hardklor parameters */
  InitStringParam("hardklor-algorithm", "version1", "basic|version1|version2",
    "Determines which spectral feature detection algorithm to use. Different results are "
    "possible with each algorithm, and there are pros and cons to each.[[html: There are "
    "three algorithms to choose from:<ul><li>basic &ndash; Performs unoptimized "
    "deconvolution and is provided for legacy purposes only.</li><li>version1 &ndash; "
    "Uses the optimizations developed during the 1.0+ series. It is very accurate, but has "
    "limited sensitivity, and moderate speed improvements.</li><li>version2 &ndash; Uses "
    "the optimizations developed for version 2.0+. It is highly sensitive, but less "
    "accurate for very low abundance features, and performs exceptionally fast.</li></ul>]]",
    "Available for crux hardklor", true);
  InitStringParam("averagine-mod", "",
    "Defines alternative averagine models in the analysis that incorporate additional "
    "atoms and/or isotopic enrichments. Modifications are represented as text strings. "
    "Inclusion of additional atoms in the model is done using by entering an atomic "
    "formula, such as: PO2 or Cl. Inclusion of isotopic enrichment to the model is done by "
    "specifying the percent enrichment (as a decimal) followed by the atom being enriched "
    "and an index of the isotope. For example, 0.75H1 specifies 75% enrichment of the first "
    "heavy isotope of hydrogen. In other words, 75% deuterium enrichment. Two or more "
    "modifications can be combined into the same model, and separated by spaces: B2 0.5B1",
    "Available for crux hardklor", true);
  InitIntParam("boxcar-averaging", 0, 0, BILLION,
    "Boxcar averaging is a sliding window that averages n adjacent spectra prior to feature "
    "detection. Averaging generally improves the signal-to-noise ratio of features in the "
    "spectra, as well as improving the shape of isotopic envelopes. However, averaging will "
    "also change the observed peak intensities. Averaging with too wide a window will "
    "increase the occurrence of overlapping features and broaden the chromatographic "
    "profiles of observed features. The number specified is the total adjacent scans to be "
    "combined, centered on the scan being analyzed. Therefore, an odd number is recommended "
    "to center the boxcar window. For example, a value of 3 would produce an average of the "
    "scan of interest, plus one scan on each side. A value of 0 disables boxcar averaging.",
    "Available for crux hardklor", true);
  InitIntParam("boxcar-filter", 0, 0, BILLION,
    "This parameter is only functional when boxcar-averaging is used. The filter will "
    "remove any peaks not seen in n scans in the boxcar window. The effect is to reduce "
    "peak accumulation due to noise and reduce chromatographic broadening of peaks. Caution "
    "should be used as over-filtering can occur. The suggested number of scans to set for "
    "filtering should be equal to or less than the boxcar-averaging window size. A value of "
    "0 disables filtering.",
    "Available for crux hardklor", true);
  InitDoubleParam("boxcar-filter-ppm", 10.0, 0.0, BILLION,
    "This parameter is only functional when boxcar-filter is used. The value specifies the "
    "mass tolerance in ppm for declaring a peak the same prior to filtering across all "
    "scans in the boxcar window.",
    "Available for crux hardklor", true);
  InitBoolParam("centroided", false,
    "Indicates whether the data contain profile or centroided peaks.",
    "Available for crux hardklor", true);
  InitStringParam("cdm", "Q", "B|F|P|Q|S",
    "Choose the charge state determination method.[[html: There are five methods to "
    "choose from:<ul><li>B &ndash; Basic method, assume all charge states are possible."
    "</li><li>F &ndash; Fast Fourier transform.</li><li>P &ndash; Patterson algorithm.</li>"
    "<li>Q &ndash; QuickCharge method, uses inverse peak distances.</li><li>S &ndash; "
    "Senko method, or combined Fast Fourier Transform and Patterson algorithm.</li></ul>]]",
    "Available for crux hardklor", true);
  InitIntParam("min-charge", 1, 1, BILLION,
    "Specifies the minimum charge state to allow when finding spectral features. It is "
    "best to set this value to the lowest assumed charge state to be present. If set higher "
    "than actual charge states that are present, those features will not be identified or "
    "incorrectly assigned a different charge state and mass.",
    "Available for crux hardklor", true);
  InitIntParam("max-charge", 5, 1, BILLION,
    "Specifies the maximum charge state to allow when finding spectral features. It is "
    "best to set this value to a practical number (i.e. do not set it to 20 when doing a "
    "tryptic shotgun analysis). If set higher than actual charge states that are present, "
    "the algorithm will perform significantly slower without any improvement in results.",
    "Available for crux hardklor", true);
  InitDoubleParam("corr", 0.85, 0, 1.0,
    "Sets the correlation threshold (cosine similarity) for accepting each predicted "
    "feature.",
    "Available for crux hardklor", true);
  InitIntParam("depth", 3, 1, BILLION,
    "Sets the depth of combinatorial analysis. For a given set of peaks in a spectrum, "
    "search for up to this number of combined peptides that explain the observed peaks. "
    "The analysis stops before depth is reached if the current number of deconvolved "
    "features explains the observed peaks with a correlation score above the threshold "
    "defined with the correlation parameter.",
    "Available for crux hardklor", true);
  InitBoolParam("distribution-area", false,
    "When reporting each feature, report abundance as the sum of all isotope peaks. The "
    "value reported is the estimate of the correct peak heights based on the averagine "
    "model scaled to the observed peak heights.",
    "Available for crux hardklor", true);
  InitStringParam("hardklor-data-file", "",
    "Specifies an ASCII text file that defines symbols for the periodic table.",
    "Available for crux hardklor", true);
  InitStringParam("instrument", "fticr", "fticr|orbitrap|tof|qit",
    "Indicates the type of instrument used to collect data. This parameter, combined with "
    "the resolution parameter, define how spectra will be centroided (if you provide "
    "profile spectra) and the accuracy when aligning observed peaks to the models.",
    "Available for crux hardklor", true);
  InitStringParam("isotope-data-file", "",
    "Specifies an ASCII text file that can be read to override the natural isotope "
    "abundances for all elements.",
    "Available for crux hardklor", true);
  InitIntParam("max-features", 10, 1, BILLION,
    "Specifies the maximum number of models to build for a set of peaks being analyzed. "
    "Regardless of the setting, the number of models will never exceed the number of peaks "
    "in the current set. However, as many of the low abundance peaks are noise or tail ends "
    "of distributions, defining models for them is detrimental to the analysis.",
    "Available for crux hardklor", true);
  InitIntParam("mzxml-filter", 1, 1, 2,
    "Filters the spectra prior to analysis for the requested MS/MS level. For example, if "
    "the data contain MS and MS/MS spectra, setting mzxml-filter = 1 will analyze only the "
    "MS scan events. Setting mzxml-filter = 2 will analyze only the MS/MS scan events.",
    "Available for crux hardklor", true);
  InitDoubleParam("mz-max", 0, 0, 10000,
    "Constrains the search in each spectrum to signals below this value in Thomsons. "
    "Setting to 0 disables this feature.",
    "Available for crux hardklor", true);
  InitDoubleParam("mz-min", 0, 0, 10000,
    "Constrains the search in each spectrum to signals above this value in Thomsons. "
    "Setting to 0 disables this feature.",
    "Available for crux hardklor", true);
  InitDoubleParam("mz-window", 4.0, 1.0, 20.0,
    "Only used when algorithm = version1. Defines the maximum window size in Thomsons to "
    "analyze when deconvolving peaks in a spectrum into features.",
    "Available for crux hardklor", true);
  InitDoubleParam("resolution", 100000, 1, BILLION,
    "Specifies the resolution of the instrument at 400 m/z for the data being analyzed.",
    "Available for crux hardklor", true);
  InitIntParam("scan-range-max", 0, 0, BILLION,
    "Used to restrict analysis to spectra with scan numbers below this parameter value. "
    "A value of 0 disables this feature.",
    "Available for crux hardklor", true);
  InitIntParam("scan-range-min", 0, 0, BILLION,
    "Used to restrict analysis to spectra with scan numbers above this parameter value. "
    "A value of 0 disables this feature.",
    "Available for crux hardklor", true);
  InitIntParam("sensitivity", 2, 0, 3,
    "Set the sensitivity level. There are four levels: 0 (low), 1 (moderate), "
    "2 (high), and 3 (max). Increasing the sensitivity will increase computation time, "
    "but will also yield more isotope distributions.",
    "Available for crux hardklor", true);
  InitDoubleParam("signal-to-noise", 1.0, 0.0, BILLION,
    "Filters spectra to remove peaks below this signal-to-noise ratio prior to finding "
    "features.",
    "Available for crux hardklor", true);
  InitIntParam("smooth", 0, 0, 21,
    "Uses Savitzky-Golay smoothing on profile peak data prior to centroiding the spectra. "
    "This parameter is recommended for low resolution spectra only. Smoothing data causes "
    "peak depression and broadening. Only use odd numbers for the degree of smoothing (as "
    "it defines a window centered on each data point). Higher values will produce smoother "
    "peaks, but with greater depression and broadening. Setting this parameter to 0 disables "
    "smoothing.",
    "Available for crux hardklor", true);
  InitDoubleParam("sn-window", 250.0, 0.0, BILLION,
    "Set the signal-to-noise window length (in m/z). Because noise may "
    "be non-uniform across a spectrum, this value adjusts the segment size "
    "considered when calculating a signal-over-noise ratio.",
    "Available for crux hardklor", true);
  InitBoolParam("static-sn", true,
    "Applies the lowest noise threshold of any sn_window across the entire mass range for a "
    "spectrum. Setting this parameter to 0 turns off this feature, and different noise "
    "thresholds will be used for each local mass window in a spectrum.",
    "Available for crux hardklor", true);
  InitBoolParam("hardklor-xml-output", false,
    "Output XML instead of tab-delimited text.",
    "Available for crux hardklor", false);
}
