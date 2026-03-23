#include "Params.h"
#include "model/objects.h"
#include "parameter.h"

using namespace std;

void Params::InitBullseyeOptions() {
  /* bullseye parameters */
  InitArgParam("MS1 spectra",
    "The name of a file from which to parse high-resolution spectra of intact peptides. "
    "The file may be in MS1 (.ms1), binary MS1 (.bms1), compressed MS1 (.cms1), or "
    "mzXML (.mzXML) format. Bullseye will search for PPIDs in these spectra.");
  InitArgParam("MS2 spectra",
    "The name of a file from which to parse peptide fragmentation spectra. The file may "
    "be in MS2 (.ms2), binary MS2 (.bms2), compressed MS2 (.cms2) or mzXML (.mzXML) format. "
    "Bullseye will assign high-resolution precursor masses to these spectra.");
  InitStringParam("hardklor-file", "",
    "Input hardklor file into bullseye",
    "Hidden option for crux bullseye.", false);
  InitDoubleParam("max-persist", 2.0, 0, BILLION,
    "Ignore PPIDs that persist for longer than this length of time in the MS1 spectra. The "
    "unit of time is whatever unit is used in your data file (usually minutes). These PPIDs "
    "are considered contaminants.",
    "Available for crux bullseye", true);
  InitBoolParam("exact-match", false,
    "When true, require an exact match (as defined by --exact-tolerance) between the "
    "center of the precursor isolation window in the MS2 scan and the base isotopic "
    "peak of the PPID. If this option is set to false and no exact match is observed, "
    "then attempt to match using a wider m/z tolerance. This wider tolerance is calculated "
    "using the PPID's monoisotopic mass and charge (the higher the charge, the smaller "
    "the window).",
    "Available for crux bullseye", true);
  InitIntParam("gap-tolerance", 1, 0, BILLION,
    "Allowed gap size when checking for PPIDs across consecutive MS1 scans.",
    "Available for crux bullseye", true);
  InitDoubleParam("bullseye-min-mass", 600, 0, BILLION,
    "Only consider PPIDs above this minimum mass in daltons.",
    "Available for crux bullseye", true);
  InitDoubleParam("bullseye-max-mass", 8000, 1, BILLION,
    "Only consider PPIDs below this maximum mass in daltons.",
    "Available for crux bullseye", true);
  InitDoubleParam("exact-tolerance", 10.0, 0, BILLION,
    "Set the tolerance (+/-ppm) for --exact-match.",
    "Available for crux bullseye", true);
  InitDoubleParam("persist-tolerance", 10.0, 0, BILLION,
    "Set the mass tolerance (+/-ppm) for finding PPIDs in consecutive MS1 scans.",
    "Available for crux bullseye", true);
  InitIntParam("scan-tolerance", 3, 0, BILLION,
    "Total number of MS1 scans over which a PPID must be observed to be considered real. "
    "Gaps in persistence are allowed by setting --gap-tolerance.",
    "Available for crux bullseye", true);
  InitDoubleParam("retention-tolerance", 0.5, 0, BILLION,
    "Set the tolerance (+/-units) around the retention time over which a PPID can be "
    "matches to the MS2 spectrum. The unit of time is whatever unit is used in your data "
    "file (usually minutes).",
    "Available for crux bullseye", true);
  InitStringParam("spectrum-format", "", "|ms2|bms2|cms2|mgf",
    "The format to write the output spectra to. If empty, the spectra will be "
    "output in the same format as the MS2 input.",
    "Available for crux bullseye", true);
  // crux pipeline options
  InitBoolParam("bullseye", false,
    "Run the Bullseye algorithm on the given MS data, using it to assign high-resolution "
    "precursor values to the MS/MS data. If a spectrum file ends with .ms2 or .cms2, matching "
    ".ms1/.cms1 files will be used as the MS1 file. Otherwise, it is assumed that the "
    "spectrum file contains both MS1 and MS2 scans.",
    "Available for crux pipeline", true);
  InitStringParam("search-engine", "tide-search", "comet|tide-search",
    "Specify which search engine to use.",
    "Available for crux pipeline", true);
  InitStringParam("post-processor", "percolator", "percolator|assign-confidence",
    "Specify which post-processor to apply to the search results.",
    "Available for crux pipeline", true);
  // create-docs
  InitArgParam("tool-name",
    "Specifies the Crux tool to generate documentation for. If the value is "
    "'list', then a list of available tools will be given. If the value is "
    "'default-params', then a default parameter file will be given."
    "If the value is 'param-table' then a table will be printed showing "
    "which parameters are associated with which commands.");
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
  // param-medic
  InitArgParam("spectrum-file",
    "File from which to parse fragmentation spectra.");
  InitBoolParam("pm-ignore-no-charge", true,
    "When parsing spectra for measurement error estimation, ignore those without charge state information.",
    "Available for param-medic, tide-search, comet, and kojak", false);
  InitDoubleParam("pm-min-precursor-mz", 400,
    "Minimum precursor m/z value to use in measurement error estimation.",
    "Available for param-medic, tide-search, comet, and kojak", true);
  InitDoubleParam("pm-max-precursor-mz", 1800,
    "Minimum precursor m/z value to use in measurement error estimation.",
    "Available for param-medic, tide-search, comet and kojak", true);
  InitDoubleParam("pm-min-frag-mz", 150,
    "Minimum fragment m/z value to use in measurement error estimation.",
    "Available for param-medic, tide-search, comet and kojak", true);
  InitDoubleParam("pm-max-frag-mz", 1800,
    "Maximum fragment m/z value to use in measurement error estimation.",
    "Available for param-medic, tide-search, comet, and kojak", true);
  InitIntParam("pm-min-scan-frag-peaks", 40,
    "Minimum fragment peaks an MS/MS scan must contain to be used in measurement error estimation.",
    "Available for param-medic, tide-search, comet, and kojak", true);
  InitDoubleParam("pm-max-precursor-delta-ppm", 50,
    "Maximum ppm distance between precursor m/z values to consider two scans "
    "potentially generated by the same peptide for measurement error estimation.",
    "Available for param-medic, tide-search, comet, and kojak", true);
  InitStringParam("pm-charges", "0,2,3,4",
    "Precursor charge states to consider MS/MS spectra from, in measurement error estimation, "
    "provided as comma-separated values.",
    "Available for param-medic, tide-search, comet, and kojak", true);
  InitIntParam("pm-top-n-frag-peaks", 30,
    "Number of most-intense fragment peaks to consider for measurement error estimation, per MS/MS spectrum.",
    "Available for param-medic, tide-search, comet, and kojak", true);
  InitIntParam("pm-pair-top-n-frag-peaks", 5,
    "Number of fragment peaks per spectrum pair to be used in fragment error "
    "estimation.",
    "Available for param-medic, tide-search, comet, and kojak", true);
  InitIntParam("pm-min-common-frag-peaks", 20,
    "Number of the most-intense peaks that two spectra must share in order to "
    "potentially be generated by the same peptide, for measurement error estimation.",
    "Available for param-medic, tide-search, comet, and kojak", true);
  InitIntParam("pm-max-scan-separation", 1000,
    "Maximum number of scans two spectra can be separated by in order to be "
    "considered potentially generated by the same peptide, for measurement error estimation.",
    "Available for param-medic, tide-search, comet, and kojak", true);
  InitIntParam("pm-min-peak-pairs", 200,
    "Minimum number of peak pairs (for precursor or fragment) that must be "
    "successfully paired in order to attempt to estimate measurement error distribution.",
    "Available for param-medic, tide-search, comet, and kojak", true);
  // localize-modification
  InitDoubleParam("min-mod-mass", 0, 0, BILLION,
    "Ignore implied modifications where the absolute value of its mass is "
    "below this value and only score the unmodified peptide.",
    "Available for localize-modification", true);

  // Kojak Parameters
  InitStringParam("auto_ppm_tolerance_pre", "false", "false|warn|fail",
    "Automatically estimate optimal value for the <code>ppm_tolerance_pre</code> "
    "parameter from the spectra themselves. false=no estimation, warn=try to "
    "estimate but use the default value in case of failure, fail=try to estimate "
    "and quit in case of failure.",
    "Available for kojak.", true);
  InitStringParam("auto_fragment_bin_size", "false", "false|warn|fail",
    "Automatically estimate optimal value for the <code>fragment_bin_size</code> "
    "parameter from the spectra themselves. false=no estimation, warn=try to "
    "estimate but use the default value in case of failure, fail=try to estimate "
    "and quit in case of failure.",
    "Available for kojak.", true);
  InitArgParam("protein database",
    "The name of the fasta file containing the amino acid protein sequences to "
    "be searched. Kojak can generate decoy sequences internally, or they may "
    "be in this file (see the <code>decoy_filter</code> option for details). It is "
    "recommended to include the full path in the name of the file.");
  InitIntParam("threads", 0,
    "Number of threads to use when searching spectra. A value of 0 will "
    "automatically match the number of threads to the number of processing "
    "cores on the computer. Additionally, negative numbers can be used to "
    "specify threads equal to all but that number of cores.",
    "Available for kojak", true);
  InitBoolParam("export_percolator", true,
    "Exports results in Percolator text format (PIN format).",
    "Available for kojak", true);
  InitBoolParam("export_pepXML", false,
    "Exports results in pepXML format.",
    "Available for kojak", true);
  InitBoolParam("export_mzID", false,
    "Exports results in mzID format.",
    "Available for kojak", true);
  InitStringParam("percolator_version", "3",
    "Changes percolator output to the format necessary for different versions of Percolator.",
    "Available for kojak", false);
  InitDoubleParam("enrichment", 0,
    "Values between 0 and 1 to describe 18O atom percent excess (APE). For "
    "example, 0.25 equals 25 APE.",
    "Available for kojak", true);
  InitIntParam("kojak_instrument", 0,
    "Values are: 0=Orbitrap, 1=FTICR (such as Thermo LTQ-FT).",
    "Available for kojak", true);
  InitBoolParam("MS1_centroid", true,
    "Are the precursor ion (MS1) scans centroided?",
    "Available for kojak", true);
  InitBoolParam("MS2_centroid", true,
    "Are the fragment ion (MS2) scans centroided?",
    "Available for kojak", true);
  InitDoubleParam("MS1_resolution", 30000,
    "Resolution at 400 m/z, value ignored if data are already centroided.",
    "Available for kojak", true);
  InitDoubleParam("MS2_resolution", 25000,
    "Resolution at 400 m/z, value ignored if data are already centroided.",
    "Available for kojak", true);
  InitStringParam("cross_link", "nK nK 138.068074 BS3",
    "Specifies the sites of cross-linking and mass modification. Four values "
    "specify a cross-link. The first two values are one or more amino acid "
    "letters (uppercase only) that can be linked. These can be the same or "
    "different depending on whether the cross-linker is homobifunctional or "
    "heterobifunctional. Use lowercase ‘n’ or ‘c’ if the linker can bind the "
    "protein termini. The third value is the net mass value of the cross-linker"
    " when bound to the peptides. The mass can be any real number, positive or "
    "negative. The identifier is any name desired for the cross-linker. If the "
    "data contain multiple cross-linkers, provide them as a comma-separated "
    "list enclosed with quotation marks. For example, a sample "
    "cross-linked with both BS3 and EDC could be specified as \"nK nK "
    "138.068074 BS3, DE nK -18.0106 EDC\".",
    "Available for kojak", true);
  InitStringParam("mono_link", "nK 156.0786",
    "Specifies the sites of incomplete cross-linking (i.e. a mono-link) and "
    "mass modification. Two values follow this parameter, separated by spaces. "
    "The first value is one or more amino acid letters (uppercase only) that "
    "can be linked. Use lowercase ‘n’ or ‘c’ if the linker can bind the protein "
    "termini. The second value is the net mass of the incomplete cross-link "
    "reaction. The mass can be any real number, positive or negative. If "
    "multiple mono-links are possible (e.g. with a heterobifunctional "
    "cross-linker), provide them as a comma-separated list enclosed in "
    "quotation marks. For example: \"nK 156.0786, nK 155.0946\".",
    "Available for kojak", true);
  InitStringParam("fixed_modification", "C 57.02146",
    "Specifies a mass adjustment to be applied to all indicated amino acids "
    "prior to spectral analysis. Amino acids are identified by their single "
    "letter designation. N-terminal and C-terminal fixed modifications are "
    "designated by n and c, respectively. The relative mass difference, "
    "positive or negative, is listed after the amino acid, separated by a "
    "space. If multiple fixed modification masses are desired, provide them as "
    "a comma-separated list enclosed in quotation marks. For example: "
    "\"C 57.02146, nK 42.01057\".",
    "Available for kojak", true);
  InitDoubleParam("fixed_modification_protC", 0,
    "Specifies a mass adjustment to be applied to all protein C-termini prior "
    "to spectral analysis. The relative mass difference may be any non-zero "
    "number.",
    "Available for kojak", true);
  InitDoubleParam("fixed_modification_protN", 0,
    "Specifies a mass adjustment to be applied to all protein N-termini prior "
    "to spectral analysis. The relative mass difference may be any non-zero "
    "number.",
    "Available for kojak", true);
  InitStringParam("modification", "M 15.9949",
    "Specifies a dynamic mass adjustment to be applied to all indicated amino "
    "acids during spectral analysis. Peptides containing the indicated amino "
    "acids are tested with and without the dynamic modification mass. Amino "
    "acids are identified by their single letter designation. N-terminal and "
    "C-terminal dynamic peptide modifications are designated by n and c, "
    "respectively. The relative mass difference, positive or negative, is "
    "listed after the amino acid, separated by a space. If multiple dynamic "
    "modification masses are desired, including to the same amino acid, provide "
    "them as a comma-separated list enclosed with quotation marks. For "
    "example: \"M 15.9949, STY 79.966331\".",
    "Available for kojak", true);
  InitStringParam("modification_protC", "0",
    "Specifies a dynamic mass adjustment to be applied to protein C-terminal "
    "amino acids during spectral analysis. Peptides containing the protein "
    "C-terminus are tested with and without the dynamic modification mass. The "
    "relative mass difference can be any non-zero value. If multiple dynamic "
    "protein C-terminal modification masses are desired, provide them as a "
    "comma-separated list enclosed in quotation marks. For example, "
    "\"56.037448, -58.005479\".",
    "Available for kojak", true);
  InitStringParam("modification_protN", "0",
    "Specifies a dynamic mass adjustment to be applied to protein N-terminal "
    "amino acids during spectral analysis. Peptides containing the protein "
    "N-terminus are tested with and without the dynamic modification mass. The "
    "N-terminus includes both the leading and 2nd amino acid, in case of "
    "removal of the leading amino acid. The relative mass difference can be "
    "any non-zero value. If multiple dynamic protein N-terminal modification "
    "masses are desired, provide them as a comma-separated list enclosed "
    "in quotation marks. For example, \"42.01055, 0.984016\".",
    "Available for kojak", true);
  InitBoolParam("diff_mods_on_xl", false,
    "Searching for differential modifications increases search "
    "times exponentially. This increase in computation can be exacerbated "
    "when searching for differential modifications on cross-linked peptides. "
    "Such computation can be avoided if is known that the cross-linked "
    "peptides should not have differential modifications. In these cases, this "
    "setting can be turned off.",
    "Available for kojak", true);
  InitIntParam("max_mods_per_peptide", 2, 0, 100000000,
    "Indicates the maximum number of differential mass modifications allowed "
    "for a peptide sequence.",
    "Available for kojak", true);
  InitBoolParam("mono_links_on_xl", false,
    "When multiple sites of linkage are available on a peptide, it is possible "
    "for that peptide to be linked to a second peptide at one site and contain "
    "a mono-link at another site. If such instances are considered rare due to "
    "the experimental conditions, then this parameter can be disabled to improve "
    "computation time.",
    "Available for kojak", true);
  InitStringParam("kojak_enzyme", "[KR] Trypsin",
    "An enzyme string code is used to define amino acid cut sites when parsing "
    "protein sequences. Following the code, a separate label can be used to "
    "name the enzyme used. The rules for peptide parsing are similar to other "
    "database search engines such as X!Tandem: 1) cleavage amino acids are "
    "specified in square braces: [], 2) a vertical line, |, indicates N- or "
    "C-terminal to the residue, 3) exception amino acids are specified in "
    "curly braces: {}.",
    "Available for kojak", true);
  InitDoubleParam("fragment_bin_size", 0.03,
    "Determines the accuracy of the scoring algorithm with smaller bins being "
    "more strict in determining matches between theoretical and observed "
    "spectral peaks. Low-resolution spectra require larger bin sizes to "
    "accommodate errors in mass accuracy of the observed peaks. Smaller bins "
    "also require more system memory, so caution must be exercised when setting "
    "this value for high-resolution spectra. For ion trap (low-res) MS/MS "
    "spectra, the recommended values is 1.0005. For high-res MS/MS, the "
    "recommended value is 0.03.",
    "Available for kojak", true);
  InitBoolParam("ion_series_A", false,
    "Should A-series fragment ions be considered?",
    "Available for kojak", true);
  InitBoolParam("ion_series_B", true,
    "Should B-series fragment ions be considered?",
    "Available for kojak", true);
  InitBoolParam("ion_series_C", false,
    "Should C-series fragment ions be considered?",
    "Available for kojak", true);
  InitBoolParam("ion_series_X", false,
    "Should X-series fragment ions be considered?",
    "Available for kojak", true);
  InitBoolParam("ion_series_Y", true,
    "Should Y-series fragment ions be considered?",
    "Available for kojak", true);
  InitBoolParam("ion_series_Z", false,
    "Should Z-series fragment ions be considered",
    "Available for kojak", true);
  InitStringParam("decoy_filter", "DECOY_ 1",
    "This parameter requires two values. The first value is a short, "
    "case-sensitive string of characters that appears in the name of every "
    "decoy protein sequence in the database. The second value is either 0 or "
    "1, where 0 indicates that these decoy sequences are already provided in the "
    "FASTA database supplied by the user, and 1 indicates Kojak should "
    "automatically generate the decoy sequences and preface the protein names "
    "with the characters supplied in the first value. If Kojak is requested to "
    "generate decoy sequences, it will save the full complement of target "
    "and decoy sequences as a fasta file in the output directory at the end of "
    "analysis. Kojak generates decoy sequences by reversing the amino acids "
    "between enzymatic cleavage sites in the protein sequence. The sites of "
    "enzymatic cleavage are determined by the rules supplied with the "
    "<code>kojak_enzyme</code> parameter. The leading methionine in the "
    "sequences, however, remains fixed. This approach ensures that, with very "
    "few exceptions, the number, length, and mass of the decoy peptides are "
    "identical to the target peptides. ",
    "Available for kojak", true);
  InitIntParam("kojak_isotope_error", 1, 0, 3,
    "Allows the searching of neighboring isotope peak masses for poorly "
    "resolve precursors. Up to three alternative isotope peak masses will be "
    "searched in addition to the presumed precursor peak mass to correct for "
    "errors in monoisotopic precursor peak identification.",
    "Available for kojak", true);
  InitIntParam("max_miscleavages", 0, 0, 100000000,
    "Number of missed enzyme cleavages allowed. If your digestion enzyme cuts "
    "at the same amino acids involved in cross-linking, then this number must "
    "be greater than 0 to identify linked peptides. In such cases, a minimum "
    "value of 2 is required to identify loop-links.",
    "Available for kojak", true);
  InitDoubleParam("max_peptide_mass", 8000.0,
    "Maximum peptide mass allowed when parsing the protein sequence database. "
    "Peptides exceeding this mass will be ignored in the analysis.",
    "Available for kojak", true);
  InitDoubleParam("min_peptide_mass", 500.0,
    "Minimum peptide mass allowed when parsing the protein sequence database. "
    "Peptides with a lower mass will be ignored in the analysis.",
    "Available for kojak", true);
  InitIntParam("min_spectrum_peaks", 12, 0, 100000000,
    "Minimum number of MS/MS peaks required to proceed with analysis of a "
    "spectrum. If spectrum_processing is enabled, the peak count occurs after "
    "the spectrum is processed",
    "Available for kojak", true);
  InitIntParam("max_spectrum_peaks", 0, 0, 100000000,
    "Maximum number of MS/MS peaks to analyze if using the "
    "<code>spectrum_processing</code> parameter. Peaks are kept in order of "
    "intensity, starting with the most intense. Setting a value of 0 keeps all "
    "peaks.",
    "Available for kojak", true);
  InitDoubleParam("ppm_tolerance_pre", 10, 0, 100000000,
    "Tolerance used when determining which peptides to search for a given "
    "MS/MS spectrum based on its precursor ion mass. The unit is "
    "parts-per-million (PPM).",
    "Available for kojak", true);
  InitBoolParam("precursor_refinement", true,
    "Some data files may filter out precursor scans to save space prior to "
    "searching. To analyze these files, the precursor analysis algorithms in "
    "Kojak must be disabled. It is also possible, though not always "
    "recommended, to disable these algorithms even when precursor scans are "
    "included in the data files. This parameter toggles the precursor analysis "
    "algorithms.",
    "Available for kojak", true);
  InitIntParam("prefer_precursor_pred", 2, 0, 2,
    "For some data (such as Thermo Orbitrap data), the MS/MS spectra may have a "
    "precursor mass prediction already. With this parameter, the Kojak algorithm "
    "can be set to either ignore, use, or supplement the predicted precursor "
    "information. There are three options for the parameter: 0 = Ignore all "
    "precursor mass predictions and have Kojak make new predictions using its "
    "precursor processing algorithms. 1 = Use the existing precursor mass "
    "predictions and skip further processing with Kojak. 2 = Use the existing "
    "precursor mass predictions and supplement these values with additional "
    "results of the Kojak precursor processing algorithms. The recommended "
    "value is 2. "
    "Supplementing the precursor values performs the following functions. "
    "First, the monoisotopic precursor mass may be refined to one determined "
    "from a point near the apex of the extracted ion chromatogram, with "
    "potentially better mass accuracy. Second, in cases where the monoisotopic "
    "peak mass might not have been predicted correctly in the original "
    "analysis, a second monoisotopic mass is appended to the spectrum, "
    "allowing database searching to proceed checking both possibilities. "
    "Third, in cases where there is obvious chimeric signal overlap, a "
    "spectrum will be supplemented with the potential monoisotopic peak masses "
    "of all presumed precursor ions. Regardless of this parameter setting, "
    "MS/MS spectra that do not have an existing precursor mass prediction will "
    "be analyzed to identify the monoisotopic precursor mass using functions "
    "built into Kojak.",
    "Available for kojak", true);
  InitBoolParam("spectrum_processing", false,
    "The MS/MS spectrum processing function will collapse isotope "
    "distributions to the monoisotopic peak and reduce the number of peaks to "
    "analyze to the number specified with the <code>max_spectrum_peaks</code> "
    "parameter.",
    "Available for kojak", true);
  InitIntParam("top_count", 20, 1, 100000000,
    "This parameter specifies the number of top scoring peptides to store in "
    "the first pass of the Kojak analysis. A second pass follows, pairing"
    "cross-linked peptides to these top sequences to produce the final "
    "cross-linked peptide score. Setting this number too low will cause "
    "cross-linked sequences to be missed. Setting this number too high will "
    "degrade the performance of the algorithm. Optimal settings will depend on "
    "database size and the number of modifications in the search. Recommended "
    "values are between 5 and 50 (20 is probably a good start).",
    "Available for kojak", true);
  InitIntParam("truncate_prot_names", 0, 0, 100000000,
    "Exports only the specified number of characters for each protein name in "
    "the Kojak output. Otherwise, if set to 0, all characters in the protein "
    "name are exported.",
    "Available for kojak", true);
  InitIntParam("e_value_depth", 5000, 1, 100000000,
    "Specifies the minimum number of tests to be present in the histogram for "
    "e-value calculations. A larger number better resolves the histogram and "
    "improves the e-value estimation for the peptide sequences in each "
    "spectrum. However, larger numbers also increase computation time. The "
    "recommended values are between 2000 and 10000",
    "Available for kojak", true);
  InitDoubleParam("min_peptide_score", 0.1, -10.0, 10.0,
    "The minimum peptide score threshold for the first (alpha) peptide during "
    "crosslink analysis. During the first pass in the analysis, if the top "
    "scoring alpha peptides do not exceed this threshold, they will not be "
    "considered for pairing with a second (beta) peptide during the second "
    "pass of the analysis. ",
    "Available for kojak", true);

  InitBoolParam("no-analytics", false, "Don't post data to Google Analytics.", "", false);
}
