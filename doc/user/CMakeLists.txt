cmake_minimum_required(VERSION 2.8.1)
cmake_policy(VERSION 2.8.1)

# Only create webpage directory if in svn workspace
if (EXISTS /noble/www/htdocs/proj/crux)
  set(WEBDIR /noble/www/htdocs/proj/crux)
  add_custom_target(
    webpage
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/doc/user
    COMMAND cp crux.html ${WEBDIR}/index.html
    COMMAND cp -f *.html ${WEBDIR}/
    COMMAND cp -f *.png ${WEBDIR}/
    COMMAND cp -f *.gif ${WEBDIR}/
    COMMAND cp *.css ${WEBDIR}/
    COMMAND cp default.params ${WEBDIR}/
    COMMAND rm -fr ${WEBDIR}/data/
    COMMAND mkdir ${WEBDIR}/data
    COMMAND cp data/demo.ms2 ${WEBDIR}/data
    COMMAND cp data/search.target.txt ${WEBDIR}/data
    COMMAND cp data/small-yeast.fasta ${WEBDIR}/data
    COMMAND cp data/target.sqt ${WEBDIR}/data
    COMMAND mkdir -p ${WEBDIR}/download
    COMMAND cp download/index.html ${WEBDIR}/download
  )
endif (EXISTS /noble/www/htdocs/proj/crux)

set (
  crux_documentation
  barista-xml-format.html
  barista.html
  broken-link.html
  bullseye.html
  bullseye_standalone_to_crux.html
  calibrate-scores.html
  comet.html
  comet.params
  compute-q-values.html 
  contributors.html
  crux-analyze-matches.html 
  crux-create-index.html 
  crux-faq.html
  crux-get-ms2-spectrum.html 
  crux-predict-peptide-ions.html 
  crux-print-processed-spectra.html
  crux-search-for-matches.html 
  crux-search-for-xlinks.html 
  crux-sequest-search.html 
  crux.css 
  crux.html 
  customize-tutorial.html 
  d.png
  data/demo.ms2 
  data/search.target.txt 
  data/small-yeast.fasta 
  data/target.sqt
  default.params 
  demos/demo1.sh
  extract-columns.html
  extract-rows.html
  fasta-format.html 
  features.html
  generate-peptides.html 
  hardklor.html
  hardklor_standalone_to_crux.html
  install-tutorial.html 
  make-pin.html
  ms2-format.html 
  ncrr.gif
  nibib.gif 
  nigms_logo-vertical.jpg
  nigms_logo.png
  nsaf_pep_equation_a.png
  nsaf_pep_equation_b.png
  nsaf_prot_equation_a.png
  nsaf_prot_equation_b.png
  parameter-file.html 
  percolator.html 
  percolator_legacy.html
  q-ranker.html 
  quant_empai_1.gif
  release-notes.html 
  sample.fasta 
  schematic.png 
  schematic.pptx
  search-tutorial.html 
  sin_pep_equation.png
  sin_prot_equation.png
  sort-by-column.html
  spectacle-eg-input.txt
  spectacle.html
  spectral-counts-tutorial.html
  spectral-counts.html
  sqt-format.html 
  start-tutorial.html 
  stat-column.html
  tide-index.html
  tide-search.html
  txt-format.html 
  UWlogo150p.gif 
)
set(
  comet_documentation
  comet-options/activation_method.html
  comet-options/add_A_alanine.html
  comet-options/add_B_user_amino_acid.html
  comet-options/add_C_cysteine.html
  comet-options/add_Cterm_peptide.html
  comet-options/add_Cterm_protein.html
  comet-options/add_D_aspartic_acid.html
  comet-options/add_E_glutamic_acid.html
  comet-options/add_F_phenylalanine.html
  comet-options/add_G_glycine.html
  comet-options/add_H_histidine.html
  comet-options/add_I_isoleucine.html
  comet-options/add_J_user_amino_acid.html
  comet-options/add_K_lysine.html
  comet-options/add_L_leucine.html
  comet-options/add_M_methionine.html
  comet-options/add_N_asparagine.html
  comet-options/add_Nterm_peptide.html
  comet-options/add_Nterm_protein.html
  comet-options/add_O_ornithine.html
  comet-options/add_P_proline.html
  comet-options/add_Q_glutamine.html
  comet-options/add_R_arginine.html
  comet-options/add_S_serine.html
  comet-options/add_T_threonine.html
  comet-options/add_U_user_amino_acid.html
  comet-options/add_V_valine.html
  comet-options/add_W_tryptophan.html
  comet-options/add_X_user_amino_acid.html
  comet-options/add_Y_tyrosine.html
  comet-options/add_Z_user_amino_acid.html
  comet-options/allowed_missed_cleavage.html
  comet-options/clear_mz_range.html
  comet-options/clip_nterm_methionine.html
  comet-options/combined.html
  comet-options/comet.params
  comet-options/content.html
  comet-options/convertCometPHPs.py
  comet-options/database_name.html
  comet-options/decoy_prefix.html
  comet-options/decoy_search.html
  comet-options/digest_mass_range.html
  comet-options/footer.html
  comet-options/fragment_bin_offset.html
  comet-options/fragment_bin_tol.html
  comet-options/head.html
  comet-options/index.html
  comet-options/isotope_error.html
  comet-options/mass_type_fragment.html
  comet-options/mass_type_parent.html
  comet-options/max_fragment_charge.html
  comet-options/max_precursor_charge.html
  comet-options/max_variable_mods_in_peptide.html
  comet-options/minimum_intensity.html
  comet-options/minimum_peaks.html
  comet-options/ms_level.html
  comet-options/nucleotide_reading_frame.html
  comet-options/num_enzyme_termini.html
  comet-options/num_output_lines.html
  comet-options/num_results.html
  comet-options/num_threads.html
  comet-options/output_outfiles.html
  comet-options/output_pepxmlfile.html
  comet-options/output_pinxmlfile.html
  comet-options/output_sqtfile.html
  comet-options/output_sqtstream.html
  comet-options/output_suffix.html
  comet-options/output_txtfile.html
  comet-options/peptide_mass_tolerance.html
  comet-options/peptide_mass_units.html
  comet-options/precursor_charge.html
  comet-options/precursor_tolerance_type.html
  comet-options/print_expect_score.html
  comet-options/remove_precursor_peak.html
  comet-options/remove_precursor_tolerance.html
  comet-options/runall
  comet-options/sample_enzyme_number.html
  comet-options/scan_range.html
  comet-options/search_enzyme_number.html
  comet-options/show_fragment_ions.html
  comet-options/skip_researching.html
  comet-options/spectrum_batch_size.html
  comet-options/theoretical_fragment_ions.html
  comet-options/topmenu.html
  comet-options/use_A_ions.html
  comet-options/use_B_ions.html
  comet-options/use_C_ions.html
  comet-options/use_NL_ions.html
  comet-options/use_X_ions.html
  comet-options/use_Y_ions.html
  comet-options/use_Z_ions.html
  comet-options/use_sparse_matrix.html
  comet-options/variable_C_terminus.html
  comet-options/variable_C_terminus_distance.html
  comet-options/variable_N_terminus.html
  comet-options/variable_N_terminus_distance.html
  comet-options/variable_mod1.html
  comet-options/variable_mod2.html
  comet-options/variable_mod3.html
  comet-options/variable_mod4.html
  comet-options/variable_mod5.html
  comet-options/variable_mod6.html
)
install (FILES ${crux_documentation} DESTINATION doc)
install (FILES ${comet_documentation} DESTINATION doc/comet-options)
