# Roadmap

There is no `TODO.md` in the repo. This roadmap uses the **95 open issues
in [crux-toolkit/crux-toolkit](https://github.com/crux-toolkit/crux-toolkit/issues)**
(snapshot 2026-10-30) as the TODO list. It orders them by scientific risk,
following `mission.md`, and splits them into very small phases. Each phase
is roughly one PR.

**Definition of done for every phase** (in addition to the phase's own
criteria):

- [ ] A failing test reproduces the issue before the fix (or the issue is
      shown to be already fixed on `master` and closed with that evidence).
- [ ] The fix is checked against an independent reference (standalone tool,
      hand calculation or known-answer data).
- [ ] Golden outputs pass on all CI platforms (once Phase 3 exists).
- [ ] `doc/user/release-notes.html` updated, including any change to output
      values.
- [ ] The PR targets `develop`, the integration branch, which is promoted to
      `master` in reviewed batches. It is stacked on another PR only when it
      builds on that PR's changes.
- [ ] The issue is closed with a link to the PR.

Legend: 🔴 silent wrong results · 🟠 crash/data loss · 🟡 build/platform ·
🔵 usability/docs · 🟢 enhancement

---

## Stage A — Safety net (no issue fixes until this exists)

### Phase 1 — Triage snapshot
- Re-run each open issue's reproduction on current `master`. Tag each one
  `confirmed`, `cannot-reproduce` or `already-fixed`.
- Likely already done: #439 (`spectrum-convert` exists), #452 (`crux-lfq`
  exists). Verify and close them.
- **Done when:** every issue has a status and a link to its repro input.

### Phase 1a — Smoke tests check exit codes · 🟠 #754
- `CruxTester#exec` always returns 0, so "the return value should be 0"
  can never fail. Return the real exit status, then fix or record every
  scenario that had been failing unnoticed.
- **Done when:** a scenario whose command fails turns red. Every existing
  scenario's exit-code expectation is verified, and those that fail are
  linked to issues.

### Phase 1b — Builds outside a git checkout don't crash · 🟠 #755
- With no `.git/HEAD`, `CRUX_VERSION` has no date and tide-search throws.
  This hits source-package builds and `git worktree` checkouts, which
  parallel phase work relies on.
- **Done when:** tide-search runs from a worktree build and from the source
  tarball, with a test covering a version string that has no date.

### Phase 2 — Unit tests build and run
- Wire `test/unit-tests` into the root CMake behind `CRUX_BUILD_TESTS`.
  Replace the PXD005590 manual download with a small fixture.
- **Done when:** `ctest` runs the crux-lfq tests locally with no manual
  steps.

### Phase 2a — Test mzML files Comet can read · 🟡 #756
- The K562 and sliced mzML files are committed with CRLF line endings, and
  Comet can't parse them. That leaves the Comet path untestable.
- **Done when:** `crux comet` runs on every committed mzML, and
  `.gitattributes` keeps them byte-exact.

### Phase 3 — Golden-output harness skeleton
- Add a pytest harness and `TOLERANCES.md` (see `tech-stack.md`). Add a
  first golden test for `tide-index` + `tide-search` on the existing
  smoke-test data.
- **Done when:** the harness passes locally and fails when a value is
  deliberately perturbed.

### Phase 4 — Scientific invariant checks
- Add generic checks to the harness: no NaN/Inf, target/decoy feature
  symmetry, monotone q-values, mass = residues + H₂O + mods. Also add the
  command-to-invariant applicability matrix (`test/golden/INVARIANTS.md`).
- **Done when:** each check runs on every golden output the matrix says it
  applies to. Any current violations are recorded as expected failures,
  each linked to its issue.

### Phase 5 — PR CI on all platforms
- A PR workflow runs build + `ctest` + golden tests on Ubuntu, macOS
  (arm64 and x86_64) and Windows.
- **Done when:** a PR is red if any platform's outputs differ beyond
  tolerance.

### Phase 6 — CI build caching · 🟡 #736
- Cache `ext/` builds keyed by dependency pins.
- **Done when:** a PR's CI time drops substantially on a warm cache.

---

## Stage B — Silent wrong results 🔴 (highest scientific risk)

### Phase 7 — Percolator target/decoy feature asymmetry · #573
- `lnNumDSP` is computed differently for targets and decoys, which breaks
  FDR.
- **Done when:** the symmetry invariant passes, and Percolator q-values
  match a hand-built pin file.

### Phase 8 — NaN retention times · #716
- **Done when:** RT is carried through correctly from the spectrum file, or
  the run fails loudly if RT is missing.

### Phase 8a — Explicit "RT unknown" instead of rt = 0 · #752
- **Done when:** every writer uses one documented missing-RT value, and every
  reader treats it as missing.

### Phase 8b — One retention-time unit across readers · #753
- **Done when:** pwiz, MSToolkit, pepXML and SQT readers all yield RT in one
  documented unit, checked on files whose RT is in seconds and in minutes.

### Phase 9 — TMT static-mod mass offset · #435, #453
- Spectrum neutral mass and peptide mass are off by multiples of the static
  mod (likely the same root cause in both issues).
- **Done when:** masses equal hand-computed values for TMT and non-TMT
  Comet/Tide inputs.

### Phase 10 — Percolator result depends on input format · #368
- **Done when:** tsv, pin, pepXML and mzIdentML inputs of the same search
  give identical Percolator output (format-equivalence invariant).

### Phase 11 — Percolator N-terminal modification handling · #375
- **Done when:** N-term mods round-trip identically through every format.

### Phase 12 — `-nan` modification mass in pepXML · #471
- **Done when:** the pepXML mod masses match the tsv output, with no NaN.

### Phase 13 — Tide ties · #388
- **Done when:** tied PSMs get equal ranks, with a deterministic, documented
  order that's identical across platforms and thread counts.

### Phase 14 — Exact p-value with non-tsv outputs · #225, #298
- **Done when:** exact-p-value results are identical across all output
  formats, and every format is available.

### Phase 15 — Percolator with exact p-values · #209
- **Done when:** the cause is explained (feature scaling or sign
  convention), fixed or documented, and benchmarked on known-answer data.

### Phase 16 — Comet adds Cys static mod implicitly · #485
- **Done when:** `crux comet` defaults match standalone Comet, or the
  difference is documented and loudly logged.

### Phase 17 — Comet variable-mod format mismatch · #731
- **Done when:** `crux comet` and standalone Comet parse the same
  `variable_modNN` strings to the same mods.

### Phase 18 — Comet decoys in target file with `concat=T` · #534
- **Done when:** target/decoy labelling is correct in every output
  configuration.

### Phase 19 — make-pin protein delimiter · #707
- **Done when:** the default delimiter matches Percolator's expectation, and
  multi-protein PSMs round-trip.

### Phase 20 — make-pin missing fields · #396
- **Done when:** every pin column is filled for every supported input
  format.

### Phase 21 — Non-standard Percolator feature files · #378
- **Done when:** extra or reordered columns are handled by name. Malformed
  files fail loudly.

### Phase 22 — pepXML splits search hits · #408
- **Done when:** all hits for one spectrum query are under a single
  `spectrum_query`, validated against the pepXML schema.

### Phase 23 — `top-match` ignored for delimited files · #303
- **Done when:** the row count per spectrum equals `top-match` in all
  formats.

### Phase 24 — Cascade-search column bug · #358
- **Done when:** the column values match a non-cascade search for the same
  PSMs.

### Phase 25 — Deisotope boundary condition · #478
- **Done when:** the `<=` vs `<` choice is justified against the algorithm
  spec, with a boundary test.

### Phase 26 — Bullseye spec / Hardklör consistency · #177, #94
- **Done when:** `crux bullseye` and `crux hardklor` match the standalone
  tools on the smoke-test data.

### Phase 27 — Protein object construction in parser · #159
- **Done when:** `MatchCollectionParser` builds the correct protein type,
  and downstream protein-level outputs are unchanged or explained.

### Phase 28 — DIAmeter duplicate fragment p-value · #645
- **Done when:** one fragment p-value remains, with its definition
  documented and golden outputs updated and explained.

### Phase 29 — Ambiguous amino-acid masses · #4
- **Done when:** B/Z/J/X handling is defined, consistent across Tide and
  Comet paths, and documented.

---

## Stage C — Crashes, data loss, wrong files 🟠

### Phase 30 — spectral-counts rejects protein database · #638
### Phase 31 — spectral-counts temp files and Windows · #409, #488
### Phase 32 — Unique temporary filenames · #491
- **Done when:** concurrent runs in the same directory never collide.

### Phase 33 — subtract-index broken · #712
### Phase 34 — Percolator FATAL writing pepXML · #469
### Phase 35 — mzIdentML read/write · #423, #429, #457
- **Done when:** output validates against the mzIdentML 1.1/1.2 schema, and
  the files in the issues parse.

### Phase 36 — cascade-search multithreading · #432
- **Done when:** results with `num-threads>1` are identical to
  single-threaded runs.

### Phase 37 — param-medic in cascade-search · #418
### Phase 38 — psm-convert crash on macOS · #318
### Phase 39 — `overwrite=T` on Windows · #442
### Phase 40 — `crux comet` with indexed databases · #533
### Phase 41 — DIAmeter output hygiene · #652, #592
- **Done when:** a run starts by clearing its own outputs, and produces only
  the documented files.

---

## Stage D — Build and platform 🟡

Stage D fixes specific breakages only. It is not a repo modernization (see
`mission.md`, Out of scope).

### Phase 42 — Fix the C++ standard macro
- `use_cxx14` sets C++11 (found during audit; no issue yet; open one). Make
  it match the C++14 actually in use. No standard upgrade.
### Phase 43 — Missing includes / removed APIs · #417, #398
### Phase 44 — Drop `-fpermissive` from Comet patch · #735
### Phase 45 — Apple Silicon build · #651
### Phase 46 — Smoke tests on macOS · #416
### Phase 47 — Lint debt · #322, then #325
### Phase 48 — Installers for macOS and Windows · #219

---

## Stage E — Test coverage gaps

### Phase 49 — Smoke tests for cascade-search · #353
### Phase 50 — Smoke tests for spectral-counts · #501
### Phase 51 — Port Ruby/Cucumber smoke tests to pytest · #670
- Split into one sub-phase per `.feature` file. Delete Cucumber and Perl at
  the end.

---

## Stage F — Provenance and usability 🔵

### Phase 52 — Version on startup and in every output · #233, #293
- **Done when:** every command logs its version, and every output format
  embeds the version and parameters (mission invariant 8).

### Phase 53 — Better "No search paths found!" error · #202
### Phase 54 — Comet index `.idx` filename · #742
### Phase 55 — Decouple `auto_fragment_bin_tol` from Comet · #705
### Phase 56 — Group parameter-file options by command · #223
### Phase 57 — assign-confidence: missing columns and generic TSV · #335, #282

---

## Stage G — Documentation and housekeeping 🔵

### Phase 58 — Home page entries · #722, #711
### Phase 59 — Tutorials: pipeline, DeepRT · #609, #593
### Phase 60 — Developer docs on usage tracking; GA4 · #346, #566
### Phase 61 — Update Crux APIs docs · #162

---

## Stage H — Enhancements 🟢 (only after Stages A–G)

Each enhancement gets golden tests before it merges, and must not change
existing outputs unless the user turns it on.

### Phase 62 — Tide index/search ion-type flags · #311, #723
### Phase 63 — Tide fragment charges above 2+ · #64
### Phase 64 — ECD/ETD spectra in Tide · #165
### Phase 65 — Selenocysteine and stop codons in tide-index · #333, #272
### Phase 66 — PEFF support in tide-index · #332, #683
### Phase 67 — Digestion rules instead of enzyme names · #110
### Phase 68 — Modification affects digestion/cross-linking flag · #85
### Phase 69 — Amino-acid-dependent neutral losses · #91
### Phase 70 — predict-peptide-ions with modifications · #389
### Phase 71 — b/y ion fractions in txt output · #672
### Phase 72 — Comet E-value as Percolator feature · #547
### Phase 73 — Q-Ranker spec features with modified peptides · #167
### Phase 74 — PinWriter uses DelimitedFileWriter · #113
### Phase 75 — param-medic in m/z or ppm · #440
### Phase 76 — spectral-counts with multiple samples · #467
### Phase 77 — DIAmeter XCorr peak reuse optimization · #548
- **Done when:** output is bit-identical to before, and it's faster.
### Phase 78 — DIAmeter multithreading · #620
- **Done when:** output is identical to single-threaded runs.
### Phase 79 — DIAmeter quantification · #665
### Phase 80 — TMT/iTRAQ quantitation command · #438
### Phase 81 — SpectrumCollection binary search · #5
- **Done when:** output is bit-identical, and it's faster.
### Phase 82 — ProteomeDiscoverer node · #164

---

## Issue coverage check

All 95 open issues appear above exactly once:
4, 5, 64, 85, 91, 94, 110, 113, 159, 162, 164, 165, 167, 177, 202, 209,
219, 223, 225, 233, 272, 282, 293, 298, 303, 311, 318, 322, 325, 332, 333,
335, 346, 353, 358, 368, 375, 378, 388, 389, 396, 398, 408, 409, 416, 417,
418, 423, 429, 432, 435, 438, 439, 440, 442, 452, 453, 457, 467, 469, 471,
478, 485, 488, 491, 501, 533, 534, 547, 548, 566, 573, 592, 593, 609, 620,
638, 645, 651, 652, 665, 670, 672, 683, 705, 707, 711, 712, 716, 722, 723,
731, 735, 736, 742.

Issues opened after the 2026-10-30 snapshot are listed here when they are
added to a phase: #754 (Phase 1a), #755 (Phase 1b), #756 (Phase 2a),
#752 (Phase 8a), #753 (Phase 8b).
