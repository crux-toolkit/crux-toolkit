#include "AminoAcidUtil.h"
#include "crux_version.h"
#include "mass.h"
#include "model/Peptide.h"
#include "model/objects.h"
#include "parameter.h"
#include "Params.h"
#include "StringUtils.h"

#include <algorithm>

using namespace std;
#include "app/CometApplication.h"
#include "app/CometIndexApplication.h"

/**
 * \file Params.cpp
 *
 * This module handles user-level parameters in Crux.  Required
 * arguments are given on the command line; optional parameters can be
 * specified either on the command line or in a parameter file.
 * Textual descriptions of each parameter are stored in the source
 * code, and these are used to automatically generate usage
 * statements, comments in output parameter files, and HTML
 * documentation.
 *
 * Following are the steps required to add a new argument or parameter
 * to an existing command:
 *
 * - In Params.cpp, inside the constructor Params::Params(), add a
 *   call to Init{Bool,Int,String,Double,Arg}Param.  This will specify
 *   things like the parameter name, type, default value, what
 *   commands it works with, and whether it is visible to the end
 *   user.
 *
 * - If your parameter is optional, then in Params::Categorize(), add
 *   an "items.insert()" with the name of your new parameter. This
 *   controls what category the parameter gets printed in when the
 *   HTML documentation is created.  Note that if the parameter
 *   doesn't get added to a category in Params::Categorize(), it just
 *   appears under a generic category called "<command name> options"
 *
 * - In the .cpp file that contains the main() for the command that
 *   uses this parameter, find a call to either getArgs() (if your new
 *   parameter is required) or getOptions() (if it is optional), and
 *   add the name of your new parameter to the list of options there.

 * - In the same file, add a call to
 *   Params::Get{Bool,Int,String,Double}() to retrieve the value
 *   associated with the parameter.  In general, these methods can be
 *   used anywhere in the source code in order to retrieve parameters.
 *   However, it's good form, when feasible, to access parameters in
 *   the main() and then to pass them as arguments, rather than
 *   accessing them as globals within subroutines.
 *
 * - If you need to edit the textual description of the command
 *   itself, search in the same file for a call to getDescription().
 *   Portions of that description enclosed in "[[nohtml: XXX ]]" will
 *   be used for the usage statement, and portions in "[[html: XXX ]]"
 *   will be in the HTML docs.
 *
 */

static Params paramContainer_;

Params::Params() : finalized_(false) {
  /* generate_peptide arguments */
  InitArgParam("protein fasta file",
    "The name of the file in FASTA format from which to retrieve proteins.");
  InitArgParam("index name",
    "The desired name of the binary index.");
  InitArgParam("ms2 file",
    "The name of one or more files from which to parse the fragmentation "
    "spectra, in any of the file formats supported by ProteoWizard.");
  /* psm-convert arguments */
  InitArgParam("input PSM file",
    "The name of a PSM file in tab-delimited text, SQT, pepXML or mzIdentML format");
  InitArgParam("output format",
    "The desired format of the output file. Legal values are tsv, html, sqt, pin, "
    "pepxml, mzidentml.");
  /* get-ms2-spectrum */
  InitArgParam("scan number",
    "Scan number identifying the spectrum.");
  InitArgParam("output file",
    "File where spectrum will be written.");
  /* predict-peptide-ions */
  InitArgParam("peptide sequence",
    "The peptide sequence.");
  InitArgParam("charge state",
    "The charge state of the peptide.");
  /* hardklor arguments */
  InitArgParam("spectra",
    "The name of a file from which to parse high-resolution spectra. The file "
    "may be in MS1 (.ms1), binary MS1 (.bms1), compressed MS1 (.cms1), or mzXML "
    "(.mzXML) format.");
  /*Percolator arguments*/
  InitArgParam("peptide-spectrum matches",
    "One or more collections of target and decoy peptide-spectrum matches (PSMs). Input may "
    "be in one of four formats: PIN, SQT, pepXML, or [[html:<a href=\"../file-formats/txt-format.html\">]]"
    "Crux tab-delimited text[[html:</a>]]. "
    "Note that if the input is provided as SQT, pepXML, or Crux "
    "tab-delimited text, then a PIN file will be generated in the output directory "
    "prior to execution. "
    "Crux determines the format of the input file by examining its "
    "filename extension."
    "[[html:<br>For PIN files, target and decoy PSMs are assumed to appear in the "
    "same file.  For other file types, decoy PSMs can be provided to Percolator in two "
    "ways: either as a separate file or embedded within the same file as the target "
    "PSMs. Percolator will first search for target PSMs in a separate file. The "
    "decoy file name is constructed from the target name by replacing \"target\" with "
    "\"decoy\". For example, if search.target.txt is provided as input, then "
    "Percolator will search for a corresponding file named search.decoy.txt. If no "
    "decoy file is found, then Percolator will assume that the given input file "
    "contains a mix of target and decoy PSMs. Within this file, decoys are identified "
    "using a prefix (specified via --decoy-prefix) on the protein name.]]");
  /*make-pin arguments*/
  InitIntParam("max-charge-feature", 0, 0, BILLION,
    "Specifies the maximum charge state feature.  When set to zero, use the "
    "maximum observed charge state.",
    "Available for make-pin and percolator.", true);
  InitArgParam("psm results",
    "A collection of target and decoy peptide-spectrum matches (PSMs). Input may be in "
    "one of four formats: SQT, PepXML (obtained from SEQUEST), [[html:<a href=\""
    "../file-formats/txt-format.html\">]]Crux tab-delimited text[[html:</a>]], or list of files (when "
    "list-of-files=T)."
    "[[html:<br>Decoy PSMs can be provided to make-pin in two ways: either as a separate "
    "file or embedded within the same file as the target PSMs. make-pin will first search "
    "for the target PSMs in a separate file. The decoy file name is constructed from the "
    "target name by replacing \"target\" with \"decoy\". For example, if search.target.txt "
    "is provided as input, then make-pin will search for a corresponding file named "
    "search.decoy.txt. If no decoy file is found, then make-pin will assume that the "
    "given input file contains a mix of target and decoy PSMs. Within this file, decoys "
    "are identified using a prefix (specified via --decoy-prefix) on the protein name.]]");
  InitStringParam("decoy input", "",
    "make-pin can convert any file format in sqt, tab-delimited and pep.xml file "
    "to pin file ",
    "Argument, not option for make-pin", false);
  InitStringParam("output-file", "",
    "Path where pin file will be written instead of make-pin.pin.",
    "It is optional for make-pin", true);
  InitBoolParam("filestem-prefixes", false,
    "Prefix PSM IDs with filestems instead of target or decoy and file index.",
    "Available for make-pin", false);
  InitBoolParam("mod-symbols", false,
    "Print modification symbols instead of masses in peptide sequences.",
    "Available for make-pin", false);
  // pipeline arguments
  InitArgParam("mass spectra",
    "The name of the file(s) from which to parse the fragmentation spectra, in any of the "
    "[[html:<a href=\"http://proteowizard.sourceforge.net/tools.shtml\">]]file formats "
    "supported by ProteoWizard[[html:</a>]]. Alteratively, with Tide-search, these files "
    "may be binary spectrum files produced by a previous run of [[html:<code>]]crux "
    "tide-search[[html:</code>]] using the [[html:<code>]]store-spectra[[html:</code>]] "
    "parameter. Multiple files can be included on the command line (space delimited), "
    "prior to the name of the database.");
  InitArgParam("peptide source",
    "Either the name of a file in fasta format from which to retrieve proteins and "
    "peptides or an index created by a previous run of [[html:<code>]]crux tide-index"
    "[[html:</code>]] (for Tide searching).");
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
  InitStringParam("require_variable_mod", "0", "0|1",
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

  // added by Yang
  /* DIAmeter-related options */
  InitStringParam("predrt-files", "",
    "The name of file from which to parse the predicted retention time of each peptide in the database. "
    "The file is tab-delimited where the first column is peptide and the second column is the predicted rt information. "
    "The rt prediction doesn't require normalization beforehand. "
    "If the peptide in the database is missing in the prediction, its predicted value will be imputed by the median of all predicted values.",
    "It is optional but recommended for DIAmeter. It can be easily generated by DeepRT or any off-the-shelf "
    "RT prediction tools by feeding in the peptide-list generated by tide-index", true);

  InitIntParam("msamanda-regional-topk", 10, 1, 1000000,
    "Analogous to the peak-picking in MS Amanda, the m/z range is divided into 10 equal length segments, "
    "and in each segment, the k most intense peaks are preserved.",
    "It is used for DIAmeter", true);

  InitIntParam("coelution-oneside-scans", 3, 1, 100,
    "(2*coelution-oneside-scans+1) scans will be used to construct the chromatogram for coelution correlation analysis",
    "It is used for DIAmeter", true);

  InitIntParam("coelution-topk", 1, 1, 10,
    "The number of topk values to consider to calculate precursor fragment co-elution",
    "It is used for DIAmeter", true);

  InitDoubleParam("coeff-precursor", 1.0, 0, 100,
    "The coefficient to balance the precursor intensity rank feature in calculating the aggregated score",
    "It is used for DIAmeter", false);

  InitDoubleParam("coeff-fragment", 1.0, 0, 100,
    "The coefficient to balance the fragment matching p-value feature in calculating the aggregated score",
    "It is used for DIAmeter", false);

  InitDoubleParam("coeff-rtdiff", 1.0, 0, 100,
    "The coefficient to balance the retention time difference feature in calculating the aggregated score",
    "It is used for DIAmeter", false);

  InitDoubleParam("coeff-elution", 1.0, 0, 100,
    "The coefficient to balance the precursor and fragment co-elution feature in calculating the aggregated score",
    "It is used for DIAmeter", false);

  InitStringParam("coeff-tag", "",
    "The tag to encode the information about coefficients",
    "It is optional but recommended for DIAmeter. If not provided, DIAmeter will automatically generate on combining all coefficients.", false);

  InitIntParam("prec-ppm", 10, 1, 1000000,
    "Tolerance used for matching precursors to spectra. "
    "Peptides must be within +/- ‘precursor-ppm’ parts-per-million (ppm) of the spectrum precursor m/z",
    "It is used for DIAmeter. Default = 10 (orbitrap) and =30 (triptof)", true);

  InitIntParam("frag-ppm", 10, 1, 1000000,
    "Tolerance used for matching fragment ions to spectrum peaks. "
    "Fragment ions must be within +/- 'fragment-ppm' of the spectrum peak value.",
    "It is used for DIAmeter. Default = 10 (orbitrap) and =30 (triptof)", true);

  InitBoolParam("unique-scannr", false,
    "Make the ScanNr of each PSM unique in the .pin file.",
    "It is used for make-pin", true);

  InitBoolParam("psm-filter", false,
    "Filter the PSM by the ensemble score.",
    "It is used for DIAmeter", true);

  InitBoolParam("spectra-denoising", false,
      "Eliminate MS2 peak if neither of the adjacent scans contains the same peak within a specified tolerance.",
      "It is used for DIAmeter", true);

  InitStringParam("diameter-instrument", "na", "orbitrap|tof5600|tof6600|na",
    "Specify the instrument platform used to acquire the input spectra. "
    "This option selects among different sets of coefficient values for the scores computed by diameter. "
    "Specifically, the 'orbitrap' setting is equivalent to spectra-denoising=false, psm-filter=false, prec-ppm=10,frag-ppm=10; "
    "the 'tof5600' is equivalent to spectra-denoising=true, psm-filter=true, prec-ppm=30, frag-ppm=30, coeff-precursor=3.2, coeff-fragment=0.2, coeff-rtdiff=0.2, coeff-elution=0.2; "
    "the 'tof6600' is equivalent to spectra-denoising=false, psm-filter=true, prec-ppm=30, frag-ppm=30, coeff-precursor=25.6, coeff-fragment=0.2, coeff-rtdiff=0.2, coeff-elution=0. "
    "Use diameter-instrument=na to set individual parameters.  Otherwise, parameters set using diameter-instrument will override any parameters set separately. ",
    "It is used for DIAmeter", true);

  Categorize();
}

Params::~Params() {
  for (map<string, Param*>::iterator i = params_.begin(); i != params_.end(); i++) {
     delete i->second;
  }
}

void Params::Categorize() {
  set<string> items;

  items.clear();
  items.insert("bullseye-max-mass");
  items.insert("bullseye-min-mass");
  items.insert("gap-tolerance");
  items.insert("max-persist");
  items.insert("persist-tolerance");
  items.insert("scan-tolerance");
  AddCategory("Identifying PPIDs in MS1 spectra", items);

  items.clear();
  items.insert("exact-match");
  items.insert("exact-tolerance");
  items.insert("retention-tolerance");
  AddCategory("Matching PPIDs to MS2 spectra", items);

  items.clear();
  items.insert("clip-nterm-methionine");
  items.insert("isotopic-mass");
  items.insert("max-length");
  items.insert("max-mass");
  items.insert("min-length");
  items.insert("min-mass");
  AddCategory("Peptide properties", items);

  items.clear();
  items.insert("cterm-peptide-mods-spec");
  items.insert("cterm-protein-mods-spec");
  items.insert("max-mods");
  items.insert("min-mods");
  items.insert("mod-precision");
  items.insert("mods-spec");
  items.insert("nterm-peptide-mods-spec");
  items.insert("nterm-protein-mods-spec");
  for (char c = 'A'; c <= 'Z'; c++) {
    items.insert(string(1, c));
  }
  items.insert("auto-modifications");
  items.insert("fixed_modification");
  items.insert("fixed_modification_protC");
  items.insert("fixed_modification_protN");
  items.insert("max_mods_per_peptide");
  items.insert("modification");
  items.insert("modification_protC");
  items.insert("modification_protN");
  AddCategory("Amino acid modifications", items);

  items.clear();
  items.insert("allow-dups");
  items.insert("decoy-format");
  items.insert("num-decoys-per-target");
  items.insert("keep-terminal-aminos");
  items.insert("seed");
  AddCategory("Decoy database generation", items);

  items.clear();
  items.insert("custom-enzyme");
  items.insert("digestion");
  items.insert("enzyme");
  items.insert("missed-cleavages");
  items.insert("kojak_enzyme");
  items.insert("max_miscleavages");
  AddCategory("Enzymatic digestion", items);

  items.clear();
  items.insert("auto-precursor-window");
  items.insert("max-precursor-charge");
  items.insert("min-precursor-charge");
  items.insert("precursor-window");
  items.insert("precursor-window-type");
  AddCategory("Precursor selection", items);

  items.clear();
  items.insert("auto-mz-bin-width");
  items.insert("compute-sp");
  items.insert("deisotope");
  items.insert("exact-p-value");
  items.insert("fragment-mass");
  items.insert("isotope-error");
  items.insert("max-ion-charge");
  items.insert("min-peaks");
  items.insert("mod-mass-format");
  items.insert("mz-bin-offset");
  items.insert("mz-bin-width");
  items.insert("peptide-centric-search");
  items.insert("remove-precursor-peak");
  items.insert("remove-precursor-tolerance");
  items.insert("scan-number");
  items.insert("skip-preprocessing");
  items.insert("spectrum-charge");
  items.insert("spectrum-max-mz");
  items.insert("spectrum-min-mz");
  items.insert("use-flanking-peaks");
  items.insert("use-neutral-loss-peaks");
  items.insert("score-function");
  items.insert("fragment-tolerance");
  items.insert("evidence-granularity");
  items.insert("top_count");
  items.insert("e_value_depth");
  items.insert("override-charges");
  AddCategory("Search parameters", items);

  items.clear();
  items.insert("picked-protein");
  items.insert("protein-enzyme");
  items.insert("protein-report-duplicates");
  items.insert("protein-report-fragments");
  AddCategory("Protein inference options", items);

  items.clear();
  items.insert("protein");
  items.insert("fido-alpha");
  items.insert("fido-beta");
  items.insert("fido-empirical-protein-q");
  items.insert("fido-fast-gridsearch");
  items.insert("fido-gamma");
  items.insert("fido-gridsearch-depth");
  items.insert("fido-gridsearch-mse-threshold");
  items.insert("fido-no-split-large-components");
  items.insert("fido-protein-truncation-threshold");
  AddCategory("Fido options", items);

  items.clear();
  // Kojak
  items.insert("cross_link");
  items.insert("mono_link");
  items.insert("mono_links_on_xl");
  items.insert("diff_mods_on_xl");
  AddCategory("Cross-linking parameters", items);

  items.clear();
  items.insert("decoy_search");
  items.insert("peff_format");
  items.insert("peff_obo");
  items.insert("decoy_filter");
  items.insert("max_peptide_mass");
  items.insert("min_peptide_mass");
  items.insert("truncate_prot_names");
  AddCategory("Database", items);

  items.clear();
  items.insert("num-threads");
  items.insert("num_threads");
  items.insert("threads");
  AddCategory("CPU threads", items);

  items.clear();
  items.insert("auto_peptide_mass_tolerance");
  items.insert("isotope_error");
  items.insert("mass_type_fragment");
  items.insert("mass_type_parent");
  items.insert("peptide_mass_tolerance");
  items.insert("peptide_mass_tolerance_lower");
  items.insert("peptide_mass_tolerance_upper");
  items.insert("peptide_mass_units");
  items.insert("precursor_tolerance_type");
  items.insert("ppm_tolerance_pre");
  items.insert("auto_ppm_tolerance_pre");
  items.insert("kojak_isotope_error");
  AddCategory("Masses", items);

  items.clear();
  items.insert("allowed_missed_cleavage");
  items.insert("num_enzyme_termini");
  items.insert("search_enzyme_number");
  items.insert("search_enzyme2_number");
  AddCategory("Search enzyme", items);

  items.clear();
  items.insert("auto_fragment_bin_tol");
  items.insert("fragment_bin_offset");
  items.insert("fragment_bin_tol");
  items.insert("theoretical_fragment_ions");
  items.insert("use_A_ions");
  items.insert("use_B_ions");
  items.insert("use_C_ions");
  items.insert("use_X_ions");
  items.insert("use_Y_ions");
  items.insert("use_Z_ions");
  items.insert("use_Z1_ions");
  items.insert("use_NL_ions");
  items.insert("auto_fragment_bin_size");
  items.insert("fragment_bin_size");
  items.insert("ion_series_A");
  items.insert("ion_series_B");
  items.insert("ion_series_C");
  items.insert("ion_series_X");
  items.insert("ion_series_Y");
  items.insert("ion_series_Z");
  AddCategory("Fragment ions", items);

  items.clear();
  items.insert("activation_method");
  items.insert("ms_level");
  items.insert("override_charge");
  items.insert("precursor_charge");
  items.insert("scan_range");
  AddCategory("mzXML/mzML parameters", items);

  items.clear();
  items.insert("clip_nterm_methionine");
  items.insert("decoy_prefix");
  items.insert("digest_mass_range");
  items.insert("equal_I_and_L");
  items.insert("mass_offsets");
  items.insert("max_duplicate_proteins");
  items.insert("max_fragment_charge");
  items.insert("max_index_runtime");
  items.insert("max_precursor_charge");
  items.insert("nucleotide_reading_frame");
  items.insert("num_results");
  items.insert("output_suffix");
  items.insert("peff_verbose_output");
  items.insert("peptide_length_range");
  items.insert("precursor_NL_ions");
  items.insert("skip_researching");
  items.insert("spectrum_batch_size");
  items.insert("text_file_extension");
  items.insert("explicit_deltacn");
  items.insert("old_mods_encoding");
  items.insert("resolve_fullpaths");
  items.insert("pinfile_protein_delimiter");
  AddCategory("Miscellaneous parameters", items);

  items.clear();
  items.insert("clear_mz_range");
  items.insert("minimum_intensity");
  items.insert("minimum_peaks");
  items.insert("remove_precursor_peak");
  items.insert("remove_precursor_tolerance");
  items.insert("kojak_instrument");
  items.insert("enrichment");
  items.insert("MS1_centroid");
  items.insert("MS2_centroid");
  items.insert("MS1_resolution");
  items.insert("MS2_resolution");
  items.insert("min_spectrum_peaks");
  items.insert("max_spectrum_peaks");
  items.insert("precursor_refinement");
  items.insert("prefer_precursor_pred");
  items.insert("spectrum_processing");
  AddCategory("Spectral processing", items);

  items.clear();
  for (int i = 1; i <= 9; i++) {
    items.insert("variable_mod0" + StringUtils::ToString(i));
  }
  for (int i = 10; i <= 15; i++) {
    items.insert("variable_mod" + StringUtils::ToString(i));
  }
  items.insert("auto_modifications");
  items.insert("max_variable_mods_in_peptide");
  items.insert("require_variable_mod");
  items.insert("protein_modlist_file");
  AddCategory("Variable modifications", items);

  items.clear();
  items.insert("add_Cterm_peptide");
  items.insert("add_Nterm_peptide");
  items.insert("add_Cterm_protein");
  items.insert("add_Nterm_protein");
  for (char c = 'A'; c <= 'Z'; c++) {
    string aaString = string(1, c);
    string aaName = AminoAcidUtil::GetName(c);
    aaName = aaName.empty() ? "user_amino_acid" : StringUtils::Replace(aaName, " ", "_");
    items.insert("add_" + aaString + "_" + aaName);
  }
  AddCategory("Static modifications", items);

  items.clear();
  items.insert("create_peptide_index");
  items.insert("create_fragment_index");
  items.insert("fragindex_max_fragmentmass");
  items.insert("fragindex_min_fragmentmass");
  items.insert("fragindex_min_ions_report");
  items.insert("fragindex_min_ions_score");
  items.insert("fragindex_num_spectrumpeaks");
  items.insert("fragindex_skipreadprecursors");
  AddCategory("Indexing", items);

  items.clear();
  items.insert("only-psms");
  items.insert("tdc");
  items.insert("search-input");
  AddCategory("General options", items);

  items.clear();
  items.insert("c-neg");
  items.insert("c-pos");
  items.insert("maxiter");
  items.insert("percolator-seed");
  items.insert("quick-validation");
  items.insert("subset-max-train");
  items.insert("test-each-iteration");
  items.insert("test-fdr");
  items.insert("train-fdr");
  items.insert("static");
  AddCategory("SVM training options", items);

  items.clear();
  items.insert("default-direction");
  items.insert("init-weights");
  items.insert("klammer");
  items.insert("output-weights");
  items.insert("override");
  items.insert("unitnorm");
  AddCategory("SVM feature input options", items);

  items.clear();
  items.insert("pm-charges");
  items.insert("pm-max-frag-mz");
  items.insert("pm-max-precursor-delta-ppm");
  items.insert("pm-max-precursor-mz");
  items.insert("pm-max-scan-separation");
  items.insert("pm-min-common-frag-peaks");
  items.insert("pm-min-frag-mz");
  items.insert("pm-min-peak-pairs");
  items.insert("pm-min-precursor-mz");
  items.insert("pm-min-scan-frag-peaks");
  items.insert("pm-pair-top-n-frag-peaks");
  items.insert("pm-top-n-frag-peaks");
  AddCategory("param-medic options", items);

  items.clear();
  items.insert("concat");
  items.insert("decoy-prefix");
  items.insert("decoy-xml-output");
  items.insert("feature-file-out");
  items.insert("file-column");
  items.insert("fileroot");
  items.insert("brief-output");
  items.insert("list-of-files");
  items.insert("mass-precision");
  items.insert("mzid-output");
  items.insert("num_output_lines");
  items.insert("output-dir");
  items.insert("output-file");
  items.insert("export_additional_pepxml_scores");
  items.insert("output_mzidentmlfile");
  items.insert("output_pepxmlfile");
  items.insert("output_percolatorfile");
  items.insert("output_sqtstream");
  items.insert("output_sqtfile");
  items.insert("output_txtfile");
  items.insert("overwrite");
  items.insert("parameter-file");
  items.insert("peptide-list");
  items.insert("pepxml-output");
  items.insert("pin-output");
  items.insert("mztab-output");
  items.insert("pout-output");
  items.insert("precision");
  items.insert("print-search-progress");
  items.insert("print_expect_score");
  items.insert("sample_enzyme_number");
  items.insert("show_fragment_ions");
  items.insert("spectrum-format");
  items.insert("spectrum-parser");
  items.insert("sqt-output");
  items.insert("store-index");
  items.insert("store-spectra");
  items.insert("temp-dir");
  items.insert("top-match");
  items.insert("txt-output");
  items.insert("use-z-line");
  items.insert("verbosity");
  items.insert("export_percolator");
  items.insert("export_pepXML");
  items.insert("export_mzID");
  AddCategory("Input and output", items);

}

bool Params::GetBool(const string& name) {
  return Require(name)->GetBool();
}

int Params::GetInt(const string& name) {
  return Require(name)->GetInt();
}

double Params::GetDouble(const string& name) {
  return Require(name)->GetDouble();
}

string Params::GetString(const string& name) {
  return Require(name)->GetString();
}

bool Params::GetBoolDefault(const string& name) {
  return Require(name)->GetBoolDefault();
}

int Params::GetIntDefault(const string& name) {
  return Require(name)->GetIntDefault();
}

double Params::GetDoubleDefault(const string& name) {
  return Require(name)->GetDoubleDefault();
}

string Params::GetStringDefault(const string& name) {
  return Require(name)->GetStringDefault();
}

const vector<string>& Params::GetStrings(const string& name) {
  Param* param = Require(name);
  if (!param->IsArgument()) {
    throw runtime_error("Parameter '" + name + "' is not an argument");
  }
  return ((ArgParam*)param)->GetStrings();
}

string Params::GetUsage(const string& name) {
  return Require(name)->GetUsage();
}

string Params::GetFileNotes(const string& name) {
  return Require(name)->GetFileNotes();
}

bool Params::IsVisible(const string& name) {
  return Require(name)->IsVisible();
}

bool Params::IsArgument(const string& name) {
  return Require(name)->IsArgument();
}

string Params::GetAcceptedValues(const string& name) {
  return Require(name)->GetAcceptedValues();
}

bool Params::IsDefault(const string& name) {
  return Require(name)->IsDefault();
}

bool Params::Exists(const string& name) {
  return paramContainer_.Get(name) != NULL;
}

void Params::Set(const string& name, bool value) {
  paramContainer_.CanModifyCheck();
  Param* param = Require(name);
  param->Set(value);
  param->ThrowIfInvalid();
}

void Params::Set(const string& name, int value) {
  paramContainer_.CanModifyCheck();
  Param* param = Require(name);
  param->Set(value);
  param->ThrowIfInvalid();
}

void Params::Set(const string& name, double value) {
  paramContainer_.CanModifyCheck();
  Param* param = Require(name);
  param->Set(value);
  param->ThrowIfInvalid();
}

void Params::Set(const string& name, const char* value) {
  Set(name, string(value));
}

void Params::Set(const string& name, const string& value) {
  paramContainer_.CanModifyCheck();
  Param* param = Require(name);
  param->Set(value);
  param->ThrowIfInvalid();
}

void Params::AddArgValue(const string& name, const string& value) {
  paramContainer_.CanModifyCheck();
  Param* param = Require(name);
  if (!param->IsArgument()) {
    throw runtime_error("Cannot add value to '" + name + "', it is not an argument");
  }
  ((ArgParam*)param)->AddValue(value);
}

void Params::Finalize() {
  paramContainer_.FinalizeParams();
}

void Params::Write(ostream* out, bool defaults) {
  if (out == NULL || !out->good()) {
    throw runtime_error("Bad file stream for writing parameter file");
  }

  *out << "# Crux parameter file (generated by Crux version " << CRUX_VERSION << ")" << endl
       << "# Full documentation available at http://cruxtoolkit.sourceforge.net/" << endl
       << "# comet_version 2016.01 rev. 1" << endl
       << "# Everything following the \'#\' symbol is treated as a comment." << endl
       << endl;

  for (vector<const Param*>::const_iterator i = Begin(); i != End(); i++) {
    string name = (*i)->GetName();
    // Print mods and Comet parameters later
    if (!(*i)->IsVisible() || name.find('_') != string::npos) {
      continue;
    }
    *out << (*i)->GetParamFileString(defaults) << endl;
  }


  // Print Comet parameters
  *out << "####################" << endl
       << "# Comet Parameters #" << endl
       << "####################" << endl;
  for (vector<const Param*>::const_iterator i = Begin(); i != End(); i++) {
    string name = (*i)->GetName();
    // Print mods and Comet parameters later
    if (!(*i)->IsVisible() || name.find('_') == string::npos) {
      continue;
    }
    *out << (*i)->GetParamFileString(defaults) << endl;
  }

  *out << "#" << endl
       << "# COMET_ENZYME_INFO _must_ be at the end of this parameters file" << endl
       << "#" << endl
       << "[COMET_ENZYME_INFO]" << endl;

  const vector<string>& cometEnzymes = get_comet_enzyme_info_lines();
  if (cometEnzymes.empty() || defaults) {
    *out << "0.  No_enzyme                      0  -          -" << endl
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
    *out << "11. Elastase                       1  ALIV       P" << endl
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
      *out << *i << endl;
    }
  }
}

map<string, Param*>::const_iterator Params::BeginAll() {
  return paramContainer_.params_.begin();
}

map<string, Param*>::const_iterator Params::EndAll() {
  return paramContainer_.params_.end();
}

vector<const Param*>::const_iterator Params::Begin() {
  return paramContainer_.paramsOrdered_.begin();
}

vector<const Param*>::const_iterator Params::End() {
  return paramContainer_.paramsOrdered_.end();
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

  // Iterate over all categories
  for (vector<ParamCategory>::const_iterator i = paramContainer_.categories_.begin();
       i != paramContainer_.categories_.end();
       i++) {
    bool any = false;
    // Iterate over each given option and check if it is in the category
    for (vector<string>::const_iterator j = options.begin(); j != options.end(); j++) {
      Param* param = Require(*j);
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
  Add(new BoolParam(name, usage, fileNotes, visible, value));
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
  Add(new IntParam(name, usage, fileNotes, visible, value, min, max));
}

void Params::InitIntParam(
  const string& name,
  int value,
  const string& usage,
  const string& fileNotes,
  bool visible
) {
  Add(new IntParam(name, usage, fileNotes, visible, value));
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
  Add(new DoubleParam(name, usage, fileNotes, visible, value, min, max));
}

void Params::InitDoubleParam(
  const string& name,
  double value,
  const string& usage,
  const string& fileNotes,
  bool visible
) {
  Add(new DoubleParam(name, usage, fileNotes, visible, value));
}

void Params::InitStringParam(
  const string& name,
  const string& value,
  const string& validValues,
  const string& usage,
  const string& fileNotes,
  bool visible
) {
  Add(new StringParam(name, usage, fileNotes, visible, value,
                      StringUtils::Split(validValues, '|')));
}

void Params::InitStringParam(
  const string& name,
  const string& value,
  const string& usage,
  const string& fileNotes,
  bool visible
) {
  Add(new StringParam(name, usage, fileNotes, visible, value));
}

void Params::InitArgParam(
  const string& name,
  const string& usage
) {
  Add(new ArgParam(name, usage));
}

Param* Params::Require(const string& name) {
  Param* param = paramContainer_.Get(name);
  if (param == NULL) {
    throw runtime_error("Parameter '" + name + "' does not exist");
  }
  return param;
}

Param* Params::Get(const string& name) {
  map<string, Param*>::iterator i = params_.find(name);
  return (i == params_.end()) ? NULL : i->second;
}

void Params::Add(Param* param) {
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

void Params::AddCategory(const string& name, const set<string>& params) {
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

void Params::FinalizeParams() {
  if (finalized_) {
    return;
  }

  for (char c = 'A'; c <= 'Z'; c++) {
    string aa = string(1, c);
    double deltaMass = GetDouble(aa);
    if (deltaMass != 0) {
      ModificationDefinition::NewStaticMod(aa, deltaMass, ANY);
    }
  }

  if (GetString("enzyme") == "no-enzyme") {
    Set("digestion", "non-specific-digest");
    Set("missed-cleavages", 500);
  }

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
  if (std::isnan(new_value)) {
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

  finalized_ = true;
}

void Params::CanModifyCheck() const {
  if (finalized_) {
    throw runtime_error("Parameters have been finalized and cannot be modified");
  }
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
string Param::GetParamFileString(bool defaultValue) const {
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
  ss << name_ << '=' << (defaultValue ? GetStringDefault() : GetString()) << endl;
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
string BoolParam::GetAcceptedValues() const { return "T|F"; }
bool BoolParam::IsDefault() const { return value_ == original_; }
bool BoolParam::GetBool() const { return value_; }
int BoolParam::GetInt() const { return IntParam::From(value_); }
double BoolParam::GetDouble() const { return DoubleParam::From(value_); }
string BoolParam::GetString() const { return StringParam::From(value_); }
bool BoolParam::GetBoolDefault() const { return original_; }
int BoolParam::GetIntDefault() const { return IntParam::From(original_); }
double BoolParam::GetDoubleDefault() const { return DoubleParam::From(original_); }
string BoolParam::GetStringDefault() const { return StringParam::From(original_); }
void BoolParam::Set(bool value) { value_ = value; }
void BoolParam::Set(int value) { value_ = From(value); }
void BoolParam::Set(double value) { value_ = From(value); }
void BoolParam::Set(const string& value) {
  try {
    value_ = From(value);
  } catch (...) {
    throw runtime_error("Invalid value for '" + name_ + "': " + "'" + value + "' "
                        "(expected boolean)");
  }
}
bool BoolParam::From(int i) { return i != 0; }
bool BoolParam::From(double d) { return d != 0; }
bool BoolParam::From(string s) {
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
string IntParam::GetAcceptedValues() const { return "<integer>"; }
bool IntParam::IsDefault() const { return value_ == original_; }
bool IntParam::GetBool() const { return BoolParam::From(value_); }
int IntParam::GetInt() const { return value_; }
double IntParam::GetDouble() const { return DoubleParam::From(value_); }
string IntParam::GetString() const { return StringParam::From(value_); }
bool IntParam::GetBoolDefault() const { return BoolParam::From(original_); }
int IntParam::GetIntDefault() const { return original_; }
double IntParam::GetDoubleDefault() const { return DoubleParam::From(original_); }
string IntParam::GetStringDefault() const { return StringParam::From(original_); }
void IntParam::Set(bool value) { value_ = From(value); }
void IntParam::Set(int value) { value_ = value; }
void IntParam::Set(double value) { value_ = From(value); }
void IntParam::Set(const string& value) {
  try {
    value_ = From(value);
  } catch (...) {
    throw runtime_error("Invalid value for '" + name_ + "': " + "'" + value + "' "
                        "(expected int)");
  }
}
int IntParam::From(bool b) { return b ? 1 : 0; }
int IntParam::From(double d) { return (int)d; }
int IntParam::From(const string& s) { return StringUtils::FromString<int>(s); }
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
    value_(value), original_(value), min_(min), max_(max) {}
void DoubleParam::ThrowIfInvalid() const {
  if (value_ < min_ || value_ > max_) {
    throw runtime_error("Value of '" + name_ + "' must be between " +
                        StringUtils::ToString(min_) + " and " +
                        StringUtils::ToString(max_));
  }
}
string DoubleParam::GetAcceptedValues() const { return "<float>"; }
bool DoubleParam::IsDefault() const { return value_ == original_; }
bool DoubleParam::GetBool() const { return BoolParam::From(value_); }
int DoubleParam::GetInt() const { return IntParam::From(value_); }
double DoubleParam::GetDouble() const { return value_; }
string DoubleParam::GetString() const { return StringParam::From(value_); }
bool DoubleParam::GetBoolDefault() const { return BoolParam::From(original_); }
int DoubleParam::GetIntDefault() const { return IntParam::From(original_); }
double DoubleParam::GetDoubleDefault() const { return original_; }
string DoubleParam::GetStringDefault() const { return StringParam::From(original_); }
void DoubleParam::Set(bool value) { value_ = From(value); }
void DoubleParam::Set(int value) { value_ = From(value); }
void DoubleParam::Set(double value) { value_ = value; }
void DoubleParam::Set(const string& value) {
  try {
    value_ = From(value);
  } catch (...) {
    throw runtime_error("Invalid value for '" + name_ + "': " + "'" + value + "' "
                        "(expected float)");
  }
}
double DoubleParam::From(bool b) { return b ? 1 : 0; }
double DoubleParam::From(int i) { return (double)i; }
double DoubleParam::From(const string& s) { return StringUtils::FromString<double>(s); }
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
    original_(value), validValues_(validValues) {
  Set(value);
}
void StringParam::ThrowIfInvalid() const {
  if (!validValues_.empty() &&
      find(validValues_.begin(), validValues_.end(), value_) == validValues_.end()) {
    throw runtime_error("Invalid value for '" + name_ + "'; must be one of <" +
                        StringUtils::Join(validValues_, '|') + ">");
  }
}
string StringParam::GetAcceptedValues() const {
  return validValues_.empty() ? "<string>" : StringUtils::Join(validValues_, '|');
}
bool StringParam::IsDefault() const { return value_ == original_; }
bool StringParam::GetBool() const { return BoolParam::From(value_); }
int StringParam::GetInt() const { return IntParam::From(value_); }
double StringParam::GetDouble() const { return DoubleParam::From(value_); }
string StringParam::GetString() const { return value_; }
bool StringParam::GetBoolDefault() const { return BoolParam::From(original_); }
int StringParam::GetIntDefault() const { return IntParam::From(original_); }
double StringParam::GetDoubleDefault() const { return DoubleParam::From(original_); }
string StringParam::GetStringDefault() const { return original_; }
void StringParam::Set(bool value) { value_ = From(value); }
void StringParam::Set(int value) { value_ = From(value); }
void StringParam::Set(double value) { value_ = From(value); }
void StringParam::Set(const string& value) {
  value_ = value != "__NULL_STR" ? value : "";
}
string StringParam::From(bool b) { return b ? "true" : "false"; }
string StringParam::From(int i) { return StringUtils::ToString(i); }
string StringParam::From(double d) { return StringUtils::ToString(d); }
//
// ArgParam
//
ArgParam::ArgParam(const string& name, const string& usage)
  : Param(name, usage, "", false), values_(vector<string>()) {}
string ArgParam::GetAcceptedValues() const { return "<string>"; }
bool ArgParam::IsArgument() const { return true; }
bool ArgParam::IsDefault() const { return false; }
bool ArgParam::GetBool() const { return BoolParam::From(GetString()); }
int ArgParam::GetInt() const { return StringUtils::FromString<int>(GetString()); }
double ArgParam::GetDouble() const { return StringUtils::FromString<double>(GetString()); }
string ArgParam::GetString() const {
  if (values_.empty()) {
    throw runtime_error("No value for argument '" + name_ + "'");
  }
  return values_.front();
}
bool ArgParam::GetBoolDefault() const { return false; }
int ArgParam::GetIntDefault() const { return 0; }
double ArgParam::GetDoubleDefault() const { return 0; }
string ArgParam::GetStringDefault() const { return ""; }
const vector<string>& ArgParam::GetStrings() const { return values_; }
void ArgParam::Set(bool value) { values_ = vector<string>(1, StringParam::From(value)); }
void ArgParam::Set(int value) { values_ = vector<string>(1, StringParam::From(value)); }
void ArgParam::Set(double value) { values_ = vector<string>(1, StringParam::From(value)); }
void ArgParam::Set(const string& value) { values_ = vector<string>(1, value); }
void ArgParam::AddValue(const string& value) { values_.push_back(value); }

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

