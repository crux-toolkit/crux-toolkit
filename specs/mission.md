# Mission

## Why this work exists

Crux (<https://crux.ms/>) is an open-source, cross-platform toolkit for
interpreting tandem mass spectrometry proteomics data: building peptide
indexes (`tide-index`, `comet-index`), searching spectra (`tide-search`,
`comet`, `kojak`, `diameter`), estimating confidence (`percolator`,
`assign-confidence`), and quantifying (`spectral-counts`, `crux-lfq`).

Results from Crux end up in published papers, biomarker studies and
clinical-research pipelines. A wrong number that looks plausible — a
shifted precursor mass, a NaN retention time, a feature computed
differently for targets than for decoys — doesn't crash anything. It
quietly changes which peptides and proteins people report, and how
confident they say they are. **Getting it slightly wrong is the worst
outcome, because nobody notices.**

The mission over the coming weeks: **close every open issue in
[crux-toolkit/crux-toolkit](https://github.com/crux-toolkit/crux-toolkit/issues)
(95 at the time of writing, 2026-10-30), most dangerous first. Every fix
must be shown to be scientifically correct, consistent and reproducible.**

## Target audience

| Audience | What they need from us |
|---|---|
| **Proteomics researchers** running Crux commands and publishing the results | Outputs they can trust. Identical results for identical inputs. Clear errors instead of silent garbage. |
| **Pipeline integrators** (core facilities; conda, Docker, Nextflow and Galaxy users) | Stable parameter names, file formats and column semantics. Any change to these is documented in the release notes. |
| **Upstream maintainers** (the Noble lab and crux-toolkit maintainers) | Small, single-purpose PRs that are easy to review. Each has a reproducing test, a release-notes entry and no unrelated churn. |
| **Method developers** benchmarking new scoring or FDR methods against Crux | Documented, deterministic behaviour that only changes when the change is announced and explained. |

## Scientific invariants (non-negotiable)

Every change must preserve these. A PR that violates one is not mergeable,
however many issues it closes.

1. **Reproducibility.** Same binary + same inputs + same parameters ⇒ same
   outputs, regardless of thread count, run order or temp-file location.
   Where randomness exists (decoy shuffling, Percolator CV splits), it is
   driven by an explicit, recorded seed.
2. **Cross-platform consistency.** Linux, macOS (Intel and Apple Silicon)
   and Windows produce identical results within the documented tolerance
   policy (see `tech-stack.md`). If they differ, that's a bug.
3. **Target/decoy symmetry.** Every score, feature and filter is computed
   the same way for targets and decoys. Any asymmetry breaks FDR control.
4. **No silent invalid values.** NaN, Inf, negative masses, empty protein
   IDs or truncated rows must never reach an output file without an error
   or a documented sentinel.
5. **Units and conventions are explicit.** Masses say whether they are
   neutral or m/z and mono or average, and whether static or variable mods
   are included. Tolerances say Da, ppm or m/z. Delimiters and modification
   notation are consistent across every Crux output format (tsv, pin,
   pepXML, mzIdentML, sqt).
6. **Format equivalence.** The same search written as tsv, pepXML or
   mzIdentML carries the same information. A downstream tool (e.g.
   Percolator) gives the same answer whichever format it reads.
7. **Agreement with reference implementations.** Where Crux wraps an
   upstream engine (Comet, Percolator, Kojak, Hardklör/Bullseye), the
   default `crux <tool>` behaviour matches the standalone tool, or the
   difference is deliberate and documented.
8. **Provenance.** Every output records the Crux version and the parameters
   that produced it.
9. **Fail loudly.** If an input can't be parsed with confidence, stop with a
   specific, actionable error. Never guess and carry on.

## How we work

- **Reproduce first.** Each issue starts with a failing test that shows the
  bug on a small, committed input. No reproduction, no fix. (If an issue
  can't be reproduced on current `master`, record that evidence on the
  issue and close it.)
- **Prove the fix scientifically, not just mechanically.** "It no longer
  crashes" isn't enough. Check outputs against an independent reference:
  the standalone upstream tool, a hand calculation, or a known-answer
  dataset.
- **Severity order.** Silent wrong numbers come first, then crashes and data
  loss, then build and platform problems, then usability and docs, then new
  features. `roadmap.md` encodes this order.
- **One issue (or one tight cluster with a shared root cause) per PR.**
  Every PR updates `doc/user/release-notes.html` (CI enforces this). It
  also calls out any change to output values, columns, defaults or
  parameter names.
- **Behaviour changes are announced.** If a fix changes numbers users have
  already published, the release notes say which commands, which
  conditions and in which direction.

## Out of scope (for this mission)

- **Modernizing or updating the repo as a whole.** That means no sweeping
  C++ standard upgrades, no codebase-wide style or lint passes, no
  dependency bumps, no build-system overhauls and no reorganizing of
  directories. We touch tooling only when an open issue or a scientific
  invariant above requires it, and only as far as that requires. Every
  unrelated change to old code is a chance to change results without
  meaning to.
- Rewrites or refactors that no open issue calls for.
- New features beyond those already requested in open issues.
- Changes to scoring or statistical methodology unless an issue documents
  that the current behaviour is wrong. Correcting a bug is in scope;
  "improving" a method is not.
