# B-E review fix campaign - full record (2026-11-06)

This report is the complete record of the fix campaign run against the
independent scientific audit review of 6 Nov 2026 ("Research methods
review", phases B, C, D and E - 136 severity-tagged findings). Every
finding has an entry here: what it claimed (code-level, plus a plain-language
explanation), what was done about it, and the test that shows the effect.
Deferrals are recorded honestly with their reasons; nothing is silently
dropped.

Companion documents (not restated here):

- STAGE_REVIEW_ADJUDICATION_2026-11-06.md - the point-by-point judgement of
  all 136 findings (agree / agree-with-modification / rebuttal) that governed
  which fixes were made. Three items ended "closed" through adjudication and
  are covered in section 3 below (B3a-01, C1-10) or under their landing
  package (C4-02).
- IMPROVEMENTS.md - the earlier audit-driven improvement log (entries 1-30);
  the fusion coordinate repair that closes B3a-01 is recorded there.

How to read an entry: each finding below gives the problem with file and
behaviour detail, a one-paragraph plain-language version, the fix, and the
test that proves the change. Package headers give the commit, the time it
landed, the files it touched and the test file added for it.

Status legend used in the 136-row map (section 5):

- landed A..landed G-f  - fixed inside this campaign, in that work package
- landed earlier        - already fixed before the campaign; closure evidence
                          cited per row
- no action needed      - the review entry was a positive citation
- partially landed      - part of the compound claim is fixed, residual named
- deferred + reason     - not fixed, with the reason grounded in the review's
                          own verdict and project constraints

Test counts in this report are derived by collecting the named test files at
report-writing time (`pytest --collect-only`), not copied from memory.

Verification figures (full derivation in section 6): the suite went from
954 passed / 119 skipped / 3 known env-only failures before the campaign to
1208 passed / 119 skipped / the same 3 known failures after the last package;
the campaign's thirteen dedicated test files now hold 254 tests between them.

---

## 1. The campaign, work package by work package

Twelve work packages were committed to the dev branch in one pass each, in
this order: A, B, C, D, E, F, G-a, G-b, G-c, G-d, G-e, G-f. The full pytest
suite was run after every package; every run was green apart from the three
pre-existing environment-dependent failures listed in section 6.

---

### Package A - negation-safe status vocabulary (commit 0dc5391, 2026-11-06 15:56)

Files: src/neoag/status_words.py (new), src/neoag/source_chain.py,
src/neoag/evidence_states.py; tests: tests/test_status_words_20261106.py
(9 tests). Suite after: 962 passed / 119 skipped / 3 known.

**B5-01 - "NOT_SUPPORTED" read as supported.**
Problem: status fields were tested with substring membership, so a value
literally saying NOT_SUPPORTED or NOT_PASSED satisfied checks for "SUPPORTED"
and "PASS" (src/neoag/source_chain.py `_truthy` and its positive-branch
callers). In plain terms, the pipeline asked "does this say supported?" and
answered yes to a field that says "not supported", turning honest negatives
into positive evidence.
Fix: a shared module src/neoag/status_words.py implements the review's F3
rule - NOT-prefixed first, then per-occurrence matching - as
`contains_positive()` / `starts_negated()`, with a MISSING_NEGATIONS set that
keeps NOT_AVAILABLE / NOT_APPLICABLE as unassessed rather than false. The
positive-branch sites (germline route, matched-normal ALT_SUPPORTED, SOMATIC,
protein evidence, PHASED, per-field orthogonal PASS tokens) now call it.
Test: tests/test_status_words_20261106.py pins NOT_SUPPORTED/NOT_PASSED ->
false, plain values -> true, and the absence-negation semantics (9 tests).

**B5-02 - dead negation classifier kept a second, wrong vocabulary.**
Problem: source_chain carried `_status_from_text` + `_NEGATIVE_TOKENS`, an
older text-classifier with its own token lists that no production code called
- a maintenance trap because a future edit could re-wire the known-bad
matching.
Fix: deleted (zero call sites repo-wide, verified before removal). The live
path is the one fixed under B5-01.
Test: removal is verified by the suite (no import breakage) plus grep audit;
the surviving behaviour is pinned by the B5-01 tests above.

**D1-03 - a normal sample's "NOT_ALT_SUPPORTED" hard-failed the event.**
Problem: evidence_states' hard-fail branch matched the substring
"ALT_SUPPORTED", which is contained in "NOT_ALT_SUPPORTED". A matched normal
that was explicitly checked and found negative triggered the hard-fail meant
for genuine alt support - an inversion of exactly the evidence it encodes.
Fix: the branch uses `contains_positive`, so the negated form no longer
fires; the positive form still hard-fails as it should.
Test: tests/test_status_words_20261106.py asserts NOT_ALT_SUPPORTED keeps the
event alive while ALT_SUPPORTED fails it.

**D1-04 - gate membership: 'PASS' inside 'NOT_PASSED'.**
Problem: same substring class in the gate logic of evidence_states - a gate
recorded as NOT_PASSED counted as passed.
Fix: contains_positive at the gate checks; the review's same-family sites
(filter FAIL, DETECTED/EXACT_MATCH, cross-assay PASS/CONCORDANT, RNA support,
APPM reject action, safety REJECT/PASS branches) were audited and fixed in
the same change.
Test: per-family regressions in tests/test_status_words_20261106.py.

**D1-05 - dead evidence-state classifier could contradict the live one.**
Problem: evidence_states.evidence_state + STATES + EvidenceState +
SUPPORTED_TOKENS was a whole parallel classification scheme referenced only
by its own test - the review flagged dead vocabularies as drift risk.
Fix: deleted (net -1 test, the one that only tested dead code).
Test: suite green after deletion; live states remain pinned by the D1-03/04
regressions.

---

### Package B - row-level provenance (commit 1ea2fb1, 2026-11-06 15:58)

Files: src/neoag/evidence_provenance.py, src/neoag/tools/runner.py,
adapters prime/bigmhc_im/deepimmuno writers and parsers,
src/neoag/immunogenicity_composite.py; tests:
tests/test_row_provenance_20261106.py (6 tests). Suite after: 968 / 119 / 3.

**C1-01 - fabricated predictor rows entered real evidence files.**
Problem: when PRIME did not return a peptide-HLA pair,
runner._prime_evidence_rows filled the gap with predict_pair_stub - a
deterministic pseudo-random score and rank written into the file that is
supposed to hold actual tool output, with nothing marking it as fake; the
non-strict missing-tool fallback could then register the whole stubbed file
as a real tool_run. Plainly: the model-invented numbers about how a peptide
would score were stored exactly like measurements and later averaged into
the composite score.
Fix: a generic row-level column evidence_row_status (real|stub) was added to
EVIDENCE_PROVENANCE_FIELDS; attach_provenance stamps it per row and producer
marks are sticky - a file-level "real" cannot launder a stubbed row back into
evidence. Missing pairs now become blank placeholder rows marked stub; a
fully stubbed row set forces stub provenance even when the caller passed a
real record; parsers keep an explicit source_file marker (e.g. "stub")
instead of overwriting it with the parsed path; run_immunogenicity_predictors
registers stub provenance when parsed output is all-stub.
Test: tests/test_row_provenance_20261106.py - missing pairs are blank stub
rows (no fabricated score present), file-level real cannot cover a stub row,
all-stub sets force stub registration (6 tests).

**C2-02 - the composite's stub filter only looked at file headers.**
Problem: immunogenicity_composite._is_stub_evidence checked only the
file-level evidence_status column, so row-level placeholders (and rows marked
invalid_for_production) were merged into the composite anyway.
Fix: row-level status is checked first; the file-level fallback now also
rejects invalid_for_production, not just plain "stub".
Test: same file - composite merges reject stub-marked and production-invalid
rows under a real file header.

---

### Package C - peptide extraction honesty (commit 48076f1, 2026-11-06 16:03)

Files: the variant->peptide extraction path (pvactools/VEP parsing, adapters,
profile plumbing); tests: tests/test_extract_bpackage_20261106.py (11 tests)
+ updated fixtures. Suite after: 980 / 119 / 3.

**B2-01 - the wildtype peptide was invented when VEP did not supply one.**
Problem: with no WildtypeProtein CSQ field, the extractor built
`wt = ref_aa * 42` - a fake protein made of the reference amino acid repeated
42 times - and cut the "wildtype" peptide from it. MT/WT comparison scoring
then ran against an invented sequence: self-similarity and specificity
comparisons were comparing real mutant sequence to fiction.
Fix: the fabrication is gone. Without plugin WT context the variant is
dropped and counted (variants_dropped_no_protein_context); every emitted row
carries wt_peptide_status WT_AVAILABLE/WT_UNAVAILABLE so downstream MT/WT
judgments cannot silently run against missing sequence. VEP-predicted
FrameshiftSequence rows still work as before.
Test: tests/test_extract_bpackage_20261106.py - no invented backbone appears
in output, dropped variants are counted, WT_AVAILABLE rows unchanged,
frameshift path still yields peptides (11 tests).

**B2-02 - FORMAT metrics attributed to the wrong ALT at multiallelic sites.**
Problem: VAF/AD were read at FORMAT index 0 regardless of which ALT allele
the chosen CSQ entry described - a peptide from ALT[1] could be given
ALT[0]'s depth and allele fraction.
Fix: metrics are read at the ALT index the picked CSQ entry describes.
Test: asymmetric AD fixture - the ALT[1] peptide must carry the ALT[1]
numbers.

**B2-03 - stop_gained variants emitted WT-identical windows, invisible to the funnel.**
Problem: the truncation-edge path returned the WT prefix (wt[:pos0]) as the
"mutant" window, so a pure stop_gained emitted a window identical to WT -
no novel sequence, but a non-zero emission, and the consequence-level funnel
had no accounting for these variants at all.
Fix: a pure stop_gained stays peptide-negative (no invented truncation-
junction source - the review explicitly rejected one), and is now counted in
variants_no_peptides_by_consequence so the funnel explains where it went.
Test: peptide count zero for stop_gained, counter incremented.

**B2-04 - tumor sample guessed as column 10.**
Problem: _resolve_tumor_sample_index fell back to assuming the tumor was
FORMAT column 10 on multi-sample VCFs - silent sample-role inversion risk
(the same class as B1-03).
Fix: it raises, quoting the sample columns actually found; single-sample VCF
is the tumor.
Test: multi-sample VCF without resolvable roles -> raised error naming
columns.

**B2-05 - modifier consequences marked missense rows COMPLEX.**
Problem: multi_aa complexity counted VEP modifier terms (splice_region etc.)
alongside protein-changing ones, so a plain missense annotated with a nearby
modifier term became COMPLEX and was dropped under exclude_multi_aa.
Fix: complexity is based only on protein-changing consequences; modifier
terms no longer drop rows.
Test: missense + splice_region row survives exclusion; true multi-aa rows
still marked.

**B2-06 - the FILTER funnel was invisible.**
Problem: rows failing FILTER were dropped silently under pass_only, so a run
could lose variants with no record of how many or why.
Fix: every FILTER value is counted (filter_status_counts) regardless of
pass_only; the flag is forwarded from the profile via
variant_peptide_pass_only (default unchanged = PASS-only, so no behaviour
change by default).
Test: counts present even with pass_only on; knob flips retention.

**B2-07 - dead sort constant; ties fell through to input order.**
Problem: the CSQ sort key carried a constant third element (dead code) and
tied entries had no key of their own beyond entry position.
Fix: dead constant removed; the rank is a two-slot tuple (MANE, CANONICAL)
and ties keep CSQ entry order through a stable sort - deterministic for a
fixed VCF, with the input-order dependence documented at the sort site.
Test: rank tuple has two slots; equal-key entries in fixed CSQ order pick the
first entry (no shuffle is asserted - reordering the input is accepted to
reorder the pick).

**B1-04 - the pipeline trusted that a plugin run added the WT column.**
Problem: after running pvactools with the WT-protein plugin, the code assumed
the produced CSQ declares WildtypeProtein; a silently failed plugin left the
whole B2-01 fabrication downstream.
Fix: run_vep_pvacseq_annotate now verifies the produced CSQ really declares
WildtypeProtein/FrameshiftSequence and raises otherwise; column status is
recorded in the run info and the extract summary (csq_wildtype_column etc.).
Test: annotated output without the declared field -> raise + summary fields.

**B1-05 - header-declared CSQ trusted without content.**
Problem: vcf_has_csq_annotations looked only at the header line; a VCF can
declare the CSQ field and carry no payload.
Fix: samples up to 500 records for an actual CSQ payload.
Test: header-only VCF now reports no annotations.

---

### Package D - clonality penalty inversion, self-similarity blank honesty (commit 3559e49, 2026-11-06 16:07)

Files: src/neoag/ccf_v2.py, adapters (easyfuse, pvactools_parser,
variant_peptide_adapter), src/neoag/sv/phase1.py, src/neoag/schemas.py,
profiles/*.toml (new knob); tests: tests/test_ccf_selfsim_dpackage_20261106.py
(13 tests).

**C4-02 - "unresolved" clonality was punished harder than measured subclonality.**
Problem: clonality_multiplier fell through unknown states to 0.65, below the
0.75 charged to a variant actually measured as subclonal - not knowing the
clonality cost a candidate more than being known to be subclonal (a penalty
inversion).
Fix, following the review's own disposition ("do not automatically raise
unresolved penalty to match subclonal without a utility model; report
confidence and missingness separately"): unresolved and vaf_only_unresolved
become neutral 1.0, and clonality_indeterminate gets its own
indeterminate_multiplier profile knob (default 0.75) so it is separable from
the subclonal bucket and tunable.
Test: tests/test_ccf_selfsim_dpackage_20261106.py - unresolved >= subclonal
(now equal-or-better), indeterminate knob default and override (13 tests).

**C6-01 - adapters wrote self_similarity_score "0.0" where nothing was computed.**
Problem: rows with no self-similarity scan got a fabricated "0.0" -
indistinguishable from a measured "no similarity at all", the best possible
score, feeding the safety gate's self-similarity logic with fiction.
Fix: the value stays blank and a new self_similarity_status column says
UNASSESSED; applied in every emitting adapter (easyfuse, pvactools parser,
variant-peptide adapter, SV phase1) plus the schema.
Test: blank value + UNASSESSED status across all four adapter routes;
measured zeros still survive as measurements.

---

### Package E - splice chain: coordinate QC, typed merges, evidence ladders (commit 685a3dd, 2026-11-06 16:17)

Files: spladder.py, pipeline.py, coordinates.py, easyquant.py,
splice_prefilter.py, splice2neo.py, splice_junction_adapter.py, irfinder.py,
consensus.py, evidence_chains.py, splice/causal_status.py (new); tests:
tests/test_splice_chain_epackage_20261106.py (27 tests).

**B3b-01 - headerless SplAdder rows parsed positionally; a truncated tail silently dropped the event.**
Problem: event rows were consumed by position with no shape check, so a row
with a surplus or truncated trailing count lost data invisibly - an event
could vanish from the candidate pool with no error and no ledger.
Fix: rows are structurally validated before parsing (even coordinate arity,
start < end, strictly monotonic exon blocks); a malformed or truncated row
now raises instead of being quietly dropped.
Test: malformed-arity and truncated-tail fixtures both raise (27-test file).

**B3b-02 - numeric columns merged by string order.**
Problem: consolidated tables merged with a lexicographic "best" fallback, so
psi/delta_psi/qvalue resolved by string order - "9" beating "10".
Fix: typed numeric merge semantics - _MAX_FIELDS / _MAXABS_FIELDS /
_MIN_FIELDS column registries plus a numeric-aware fallback - so maxima,
max-absolutes and minima each follow declared numeric intent.
Test: qvalue 0.9 vs "0.10" style cases resolve numerically.

**B3b-03 - ImmunoPepper only saw SplAdder junctions.**
Problem: target junctions came from SplAdder output only; junctions found
exclusively by RegTools lost their translation metadata - real novel
junctions dropped from peptide building for provenance reasons, not biology.
Fix: target junctions are assembled from every junction source;
RegTools-only junctions keep their translation metadata.
Test: RegTools-only junction reaches the ImmunoPepper input with metadata.

**B3b-04 - coordinate-system claims trusted over content; STAR aliases unmapped.**
Problem: a BED (0-based) file declared under a 1-based system was accepted
blindly, and STAR-SJ/STAR_sj aliases did not reach the SJ.out.tab column
map - the latter silently zeroing STAR read counts and strand.
Fix: a content-vs-claim guard (sniff_junction_table_format +
check_junction_format_vs_system) rejects mismatches, and the aliases now map
to the SJ.out.tab columns.
Test: BED content under a 1-based claim raises; alias table parses reads and
strand.

**B3b-05 - a negative targeted requant result produced RNA cis support.**
Problem: a targeted-requant run measuring a negative (no support) delta fed
the DNA_RNA_CIS_SUPPORTED status - absence presented as positive evidence.
Fix: a negative result yields DNA_PREDICTION_ONLY; absence stays absence.
Test: negative-delta row gets DNA_PREDICTION_ONLY, never CIS_SUPPORTED.

**B3b-06 - ORF grade O3 locked out of the pass list; INCOMPLETE QC acted as contradiction.**
Problem: two operational traps in splice_prefilter: the CREDIBLE_ORF pass
whitelist {O1, O2} omitted the best grade O3 (an inverted ordinal vocabulary
- the best ORF could sit in REVIEW forever), and FORMAL_JUNCTION_READ_QC
INCOMPLETE counted as FAIL, so STAR-only data under a complete policy
rejected every splice event as if contradicted.
Fix: O3 is in the credible-ORF gate; INCOMPLETE junction-read QC maps to
REVIEW by default with an opt-in reject policy for strict deployments, so
insufficient QC reads as "insufficient", not "contradicted".
Test: O1/O2/O3 pass, worse grades do not; STAR-only + complete policy ->
REVIEW by default, REJECT under the opt-in.

**B3b-07 - "NMD_ESCAPE" matched the NMD fail token; top-N truncated by file order.**
Problem: the substring "NMD" inside the positive status NMD_ESCAPE matched
the FAIL token list, inverting an escape-from-NMD call into an NMD failure;
and per-event truncation took the first 20 rows in file order rather than
the best binders, making candidate recall order-dependent.
Fix: precedence tokens stop NMD_ESCAPE from matching the FAIL token;
per-event top-N is binder-rank ordered.
Test: NMD_ESCAPE passes the escape check; shuffled file order yields the
same top-N set.

**B3b-09 - junction-crossing and novel-aa flags defaulted to true.**
Problem: crosses_junction and contains_novel_aa returned presumed-true when
they could not tell - "unknown" pre-charged as "yes, novel".
Fix: both default to UNASSESSED (coordinates.py / splice2neo.py).
Test: unknown-input calls return UNASSESSED.

**B3b-10 - junction provenance hardcoded "RegTools".**
Problem: splice_junction_adapter entity_rows listed primary_tools as a
hardcoded string, so events supported by other tools were credited to a tool
that may never have run.
Fix: primary_tools derived from the actual evidence index.
Test: non-RegTools evidence yields the right tool name.

**B3b-11 - IRFinder site depth could satisfy the split-read gate.**
Problem: IRFinder-S site coverage was written into the split-read fields,
letting intron-retention depth count as splice-junction split-read support -
one evidence layer substituting for another.
Fix: site coverage lands in a dedicated ir_site_coverage column; split-read
fields stay blank with a split_read_support_status marker.
Test: IR-only evidence does not pass the split-read gate.

**B3b-12 - a predicted value alone earned "favorable presentation".**
Problem: the consensus presentation tier graded rows favorable on predicted
immunogenicity values without the percentile check the doc-level claim
promises (P2).
Fix: favorable requires percentile <= 2.0 (P2); weak predicted values grade
P1 with CAP_PRESENTATION_PREDICTED_NOT_FAVORABLE_R3; missing stays P0.
Test: predicted-but-not-percentile row caps at P1 with the reason code.

**B3b-13 - three drifting copies of the causal-status ladder.**
Problem: consensus.py, evidence_chains.py and pipeline.py each kept a local
priority dict for causal status; they were already free to diverge.
Fix: one shared ladder in splice/causal_status.py imported by all three.
Test: all three modules agree on ordering cases (27 tests total in the
package file).

Side effect recorded honestly: ImmunoPepper k-mer orphans now emit one
explicit summary conflict instead of being silently suppressed, so a
sequence-identity mismatch is visible per run.

---

### Package F - SV chain and RNA allele counting (commit 6c687e9, 2026-11-06 16:32)

Files: src/neoag/sv/sv_callset.py, sv_merge.py, protein_reconstruct.py,
phase1.py, score_pipeline.py, src/neoag/reports_dual.py,
scripts/rna_allele_counts_pysam; tests: tests/test_sv_rna_fpackage_20261106.py
(22 tests) + rewritten tests/test_rna_allele_counts_script.py (8 tests).
Suite after: 1048 / 119 / 3.

**B4-01 - caller support keys unparsed; absent support shown as measured zero.**
Problem: the SV parser did not know GRIDSS' FORMAT keys (SUR/SU, PAI/PAIR,
PSC/PSR, VD) or SvABA's (SPLIT_READS/PAIRED_READS/GENOME_DEPTH), so support
reads that were present came out empty - and empty was later presented as a
measured zero, i.e. "the caller saw no support" when the real state was
"we never looked".
Fix: support parsing for both callers plus a per-sample support status
vocabulary MEASURED / NO_EVIDENCE_KEYS / FORMAT_ABSENT / SAMPLE_ABSENT, so
absence can never be rendered as a measured zero.
Test: SVRecord fixtures with each key family parse; a sample without FORMAT
gets FORMAT_ABSENT, not zero (22-test file).

**B4-06 - FILTER="." (no decision made) was silently dropped.**
Problem: a VCF record whose FILTER field is "." records that no filter
decision was made; under the default PASS-only policy those records were
dropped with no record at all - the funnel showed neither how many nor that
a "no-decision" class existed.
Fix: FILTER="." records are excluded by default but written to a ledger
(parsed/sv_vcf_filter_exclusions.tsv, written every run) with reason
EXCLUDED_FILTER_NOT_ASSESSED, and kept when unfiltered_policy=include;
genuinely-failed records stay excluded as before.
Test: dot-FILTER record excluded + ledgered by default, retained under the
include policy, failed record excluded in both.

**B4-07 - a single-sample VCF could not prove somatic, but looked like it did.**
Problem: tumor specificity treated "normal had zero support" (measured) and
"there was no normal to look at" (unassessed) identically, so single-sample
SV calls carried full somatic confidence.
Fix: tumor_specificity_status distinguishes matched-normal-zero-support from
unassessed normal; the unassessed case is capped via the profile knob
normal_unassessed_specificity_cap (default 0.85); support statuses and
breakpoint provenance are exposed on the event rows.
Test: single-sample event gets the unassessed status and the cap; matched
zero-support keeps full specificity.

**B4-08 - fusion joins ignored breakend orientation.**
Problem: protein reconstruction paired breakends without reading the VCF 4.2
ALT orientation operators, so it could join a 5' partner's tail to a 3'
partner's tail - a chimera the rearrangement cannot produce.
Fix: partner roles and required strands are parsed from the canonical ALT
syntax; orientation-impossible joins are rejected before peptide generation;
symbolic ALTs (no explicit orientation) stay ORIENTATION_UNASSESSED with low
confidence instead of pretending to be resolved.
Test: forward/reverse ALT fixtures yield the biologically valid join only;
symbolic ALT -> ORIENTATION_UNASSESSED.

**B4-09 - breakpoints snapped to the nearest CDS base recorded as exact.**
Problem: projecting a breakpoint to the nearest coding base was recorded as
an exact projection, so a junction off by hundreds of bases carried exact
confidence.
Fix: projection is recorded EXACT / APPROXIMATE_nBP / UNMAPPED, and any
approximate projection forces low reconstruction confidence.
Test: off-CDS breakpoint yields APPROXIMATE_nBP and capped confidence.

**B4-10 - SV clustering depended on input order.**
Problem: greedy clustering merged calls in file order, so shuffling the input
could change the clusters - identical data, different candidates.
Fix: clustering is now the transitive closure of pairwise compatibility
(union-find over chrom-pair buckets with an endpoint sweep); the partition
is identical for every input order.
Test: permuted inputs yield identical partitions.

**B4-14 - predictor provenance followed a global flag, not the actual producer.**
Problem: a user-supplied real MHCflurry table was labelled stub because the
run's global binding_stub flag said so, while a pipeline-generated stub file
could be discovered at runtime and upgraded to real - provenance lied in
both directions.
Fix: provenance follows the actual producer per file - supplied tables stay
real under binding_stub=True, generated stub output stays stub - and runtime
output discovery no longer upgrades an honestly-registered status.
Test: both directions pinned (supplied-real stays real, generated-stub stays
stub, discovery cannot upgrade).

**B4-16 - a gene missing from the expression table got TPM 0.**
Problem: an absent row in the expression table was encoded as event_expression
0 - "this fusion's gene is not expressed" when the truth was "this gene was
not in the table we were given".
Fix: expression_evidence_status vocabulary; event_expression stays blank when
unassessed.
Test: gene absent -> blank + status; gene present at 0 -> measured zero.

**C3-01 - RNA counting ran without BAM index or contig cross-check.**
Problem: the pysam allele-counting script started iterating and crashed deep
in (or produced credible-looking zero-depth rows) when the .bai was missing
or the VCF and BAM used different contig names / builds.
Fix: fail-fast checks - .bai presence and VCF/BAM contig-name intersection -
before iteration.
Test: missing index and cross-build contigs raise explicit errors (8-test
rewritten file).

**C3-02 - only the first ALT was counted at multiallelic sites.**
Problem: pileup counted ALT[0] only, undercounting depth attributable to
other ALTs and skewing allele fractions.
Fix: per-ALT rows sharing one depth = ref + all ALTs + other.
Test: two-ALT site: each ALT row counts its own reads against the shared
depth.

**C3-03 - per-site pileup exceptions were swallowed.**
Problem: a pileup failure on a site was caught and continued, leaving a row
with zero counts that looked measured.
Fix: errors surface as ERROR rows with blank numerics; complex-ALT sites are
labelled UNASSESSED_COMPLEX.
Test: forced pileup error -> ERROR row with blank numerics, never zeros.

---

### Package G-a - runtime config, sample roles, genome-build keys (commit 09b0d88, 2026-11-06 16:35)

Files: src/neoag/utils.py (write_toml), scripts/run_candidate_upstream.py,
src/neoag/open_neo/tool_consensus.py, src/neoag/cross_site_rna.py; tests:
tests/test_gpackage_ga_20261106.py (12 tests).

**B1-02 - runtime TOML built by string interpolation.**
Problem: the candidate-upstream runtime config was assembled by pasting paths
into a template; a path containing a quote, backslash or newline silently
changed the document's meaning - the analysis would read a different path
than the one configured, or fail obscurely.
Fix: utils.write_toml emits typed TOML from data structures with a tomllib
round-trip check, so any path that does not survive serialization fails
loudly before execution.
Test: paths with quotes/backslashes/spaces round-trip; malformed output
raises (12-test file).

**B1-03 - tumor/normal roles guessed from column order.**
Problem: an ambiguous multi-sample VCF fell through role resolution to
"assume first column is normal, second is tumor" - inverting somatic evidence
silently if the VCF was built the other way round.
Fix: resolve_vcf_sample_roles fails closed when roles cannot be resolved
(explicit roles or trusted metadata required); a single-sample VCF is the
tumor.
Test: swapped-column and ambiguous-name fixtures raise instead of guessing.

**B5-03 - GRCh38 hardcoded into consensus splice keys.**
Problem: tool_consensus stamped "GRCh38" into splice source keys regardless
of the run's declared build - hg38 and hg19 evidence could be merged under a
false label.
Fix: _write_splice takes the declared run build (inputs.genome_build) as the
default for sources without a build column; row-level build columns still
win.
Test: declared-build run stamps the declared build; row column overrides.

**C6-06 - cross-site RNA keys ignored genome build.**
Problem: cross_site_rna.event_keys identified events by coordinate only, so
the same coordinates on different builds became the same event and were
compared as if equivalent.
Fix: SNV/INDEL keys embed each row's genome_build (falling back to the
declared run build); the summary records the effective default build.
Test: same coordinates under two builds produce two distinct keys.

---

### Package G-b - predictor evidence honesty (commit 9d84f14, 2026-11-06 16:46)

Files: adapters/iedb_immunogenicity.py, adapters/mhcflurry.py,
adapters/netmhcpan.py, adapters/peptide_netmhcpan.py, src/neoag/presentation.py,
src/neoag/schemas.py, src/neoag/tools/runner.py; tests:
tests/test_gpackage_gb_20261106.py (16 tests).

**C1-02 and C1-09 - absence persisted as fake worst-case sentinels.**
Problem: where a predictor produced no value, the adapters persisted 99.0 /
0.0 / "99" / "0" sentinels - fabricated "definitely not a binder" numbers
that a later reader cannot distinguish from measurements (99.0 percentile
looks like a measured rejection).
Fix: absence is persisted as "" everywhere across the mhcflurry /
netmhcpan / presentation routes; measured zeros survive untouched.
Test: parametrized absence-vs-measured-zero matrix over all writers (16-test
file; the docstring names the item each test locks).

**C1-03 - the local IEDB-style heuristic polluted tool-named columns.**
Problem: a home-grown sequence heuristic wrote into the same evidence columns
as real tool output, so a heuristic guess was indistinguishable from a
predicted measurement by column name alone.
Fix: the heuristic is off by default (opt-in), never enters the tool-named
evidence column, and when explicitly enabled lands in separate
*_local_estimate columns.
Test: default run leaves tool columns untouched; enabled run fills only the
_local_estimate pair.

**C1-04 - "NA" placeholders for missing annotation.**
Problem: missing annotation fields were filled with the string "NA", which
downstream truthiness treats as a value.
Fix: "" (typed blank) is used; paired with the row status columns from
packages B/G-b.
Test: annotation-missing rows carry "" and the status marker.

**C1-05 - one odd residue crashed a whole batch.**
Problem: a non-standard residue (X, U, B, Z, *) hit IMMUNOSCALE[aa] as a
bare KeyError and killed the batch run.
Fix: nonstandard_residues() detects them; the affected row blanks its own
IEDB row with a reason instead of crashing; the checked call raises a clean
ValueError naming the residue.
Test: "AXU*BZ" listing, batch continues with one blanked row, clean error
text.

**C1-06 - stored tool outputs reused without checking the input.**
Problem: an existing tool output file was reused whenever present, even if
the peptide list that produced it had changed - stale predictions presented
as current.
Fix: a sidecar fingerprint over peptides + tool version + mode is written
with each stored output; reuse requires a match, otherwise re-run.
Test: same file + changed peptides -> fingerprint mismatch -> recompute.

**C1-08 - the first requested allele was fabricated as the binder.**
Problem: when no predictor result existed, the code recorded the first
requested HLA allele as the "restricting allele" - inventing presentation
evidence for an allele that was never scored.
Fix: no predictor result means the restricting HLA stays unknown.
Test: empty predictor set -> allele field blank/unknown, never the first
request.

**C1-11 - trained vs pseudoallele support recorded nowhere at predict time.**
Problem: whether NetMHCpan actually trained on an allele or used a
pseudo-allele substitute was discarded at parse time, so the presentation
columns had to say UNASSESSED even when the tool had told us.
Fix: allele support status (trained / pseudoallele / unassessed) is captured
at predict time and carried through to the presentation columns.
Test: trained and pseudoallele fixtures carry through; unassessed stays
unassessed.

---

### Package G-c - composite and APM honesty (commit ffa1e37, 2026-11-06 19:31)

Files: src/neoag/immunogenicity_composite.py, src/neoag/appm_v2.py,
src/neoag/immune_escape.py; tests: tests/test_gpackage_gc_20261106.py
(24 tests).

**C2-01 - PRIME absolute scores mixed into the composite.**
Problem: the composite used PRIME's raw score column where the documented
contract is a rank-based component, so an absolute number from one tool's
scale was averaged with differently-scaled evidence; score-only rows were
also treated as assessed.
Fix: the PRIME component is rank-only; rows with a score but no rank stay
unassessed instead of being silently used.
Test: rank present -> component used; score-only row -> unassessed (24-test
file).

**C2-03 - evidence merges joined on raw HLA strings.**
Problem: merging evidence tables keyed on (peptide, HLA) text, so "A*02:01"
in one file and "HLA-A*02:01" in another failed to join - real evidence
silently lost by notation.
Fix: canonical (peptide, HLA) keys on both sides of every evidence merge.
Test: cross-notation tables join; case differences join.

**C2-04 - composite evidence tables reused without an input fingerprint.**
Problem: existing *_evidence.tsv files were picked up whenever present,
regardless of whether the inputs that produced them changed (the C1-06
pattern one layer up).
Fix: a matching input fingerprint sidecar is required; stale evidence is
re-run or registered missing - never reused quietly.
Test: fingerprint mismatch forces recompute; matching sidecar allows reuse.

**C2-06 - string-based neutral matching and a mislabelled mean.**
Problem: neutral/sentinel detection compared strings (a numeric "0" could
behave as text) and mode="mean" did not compute an arithmetic mean.
Fix: neutral-sentinel matching is numeric; mode="mean" is the actual mean.
Test: numeric zero vs blank handled distinctly; mean values verified.

**C5-01 - every non-clean LOH state capped peptides as immune escape.**
Problem: the HLA-LOH cap fired on any non-clean state including UNASSESSED
and SINGLE_TOOL_RETAINED - a missing assessment or a lone-tool retention was
narrated as the tumour losing the allele, i.e. absence spoken of as biology.
Fix: only loss-flavoured states cap; UNASSESSED and SINGLE_TOOL_RETAINED
record a reason without penalising the peptide.
Test: per-state matrix - loss caps, unassessed/retained do not, reasons
recorded.

**C5-02 - provided-but-unmapped proteomics inputs counted as gene evidence.**
Problem: a proteomics/ligandome file supplied but whose gene identifiers did
not map was effectively counted as coverage of nothing at all.
Fix: such inputs are labelled input_only_coverage_unresolved and excluded
from gene-level evidence counts.
Test: unmapped-but-supplied file -> input_only label, gene evidence unchanged.

**C5-05 - paths advertised explainability outputs before they existed.**
Problem: the paths record advertised explainability files that had not been
written, so consumers (and audits) saw outputs that were not there.
Fix: paths advertise explainability outputs only after they exist.
Test: fresh run advertises none; after build, the written ones appear.

---

### Package G-d - deterministic CN, strict purity, FASTA access, SV evidence (commit 2c69202, 2026-11-06 19:45)

Files: src/neoag/ccf_v2.py, src/neoag/sv/{reference,phase1,evidence,peptide_builder,score_pipeline,sv_callset}.py,
adapters; tests: tests/test_gpackage_gd_20261106.py (35 tests) + one
pre-existing red integration test in tests/test_sv_score.py fixed as a side
item (it was already failing at HEAD before this package - proven via a
worktree check - because it did not supply run-matched core-predictor
evidence).

**C4-01 - overlapping CNV segments resolved by file order.**
Problem: find_cn picked whichever overlapping segment came first in the file,
so copy number (and through it CCF and clonality) depended on input ordering,
and contradictory overlapping calls passed unnoticed.
Fix: segments are ranked deterministically (confidence, then width, then
position), the chosen span is recorded, and an overlap-consistency status is
emitted - DISCORDANT demotes cnv_confidence so contradictory calls cannot
silently carry high confidence.
Test: permuted segment files pick the same segment; discordant overlap
demotes confidence (35-test file).

**C4-03 - "SVclone resolution" overclaimed.**
Problem: the field named resolution recorded only that an external tool's
output existed - advisory information presented as a resolved quantity.
Fix: renamed to external_available_advisory; the name now says what it is.
Test: column presence + name pinned by fixture output check.

**C4-04 - purity range strings parsed loosely.**
Problem: _purity_bounds accepted prose ("approximately 0.4") and
percent-scaled text, silently mis-scaling purity - the single most leveraged
CCF nuisance parameter.
Fix: structured keys parse first; free text only when the whole string is
exactly a numeric pair (x-y or tuple forms); prose and percent forms are
rejected with an explicit error.
Test: prose and "40%" inputs raise; structured and "0.3-0.5" forms parse.

**B4-02 - whole genome loaded into RAM for region fetches.**
Problem: FastaReference read entire FASTA files into memory to slice regions
- an operational memory hazard on real genomes (the review tagged it
performance; it becomes a correctness risk when the process is OOM-killed
mid-run).
Fix: byte-level random access via .fai index when present, otherwise a
streaming offset scan; no whole-genome load.
Test: region extraction correctness with and without a .fai sidecar; no
full-file read (ast-level check in the test).

**B4-03 - no ligand catalogue scored as zero overlap.**
Problem: normal_hla_ligand_overlap returned 0 when no catalogue was supplied
- "this peptide looks nothing like normal HLA ligands" from an absence of
data.
Fix: reports not_assessed when no catalogue is given.
Test: catalogue-absent call returns not_assessed, not 0.

**B4-05 - duplicate dict keys and fixed priors dressed as measurements.**
Problem: the SV raw-event builder assembled its row dict with duplicated
literal keys (later silently overriding earlier ones - a latent data-loss
bug) and stamped fixed-prior clonality/persistence numbers into measured
fields.
Fix: duplicate keys removed; a new event_component_status column honestly
labels the fixed-prior components so downstream code can tell prior from
observation (the review's B3a-03 disposition: keep priors separate from
measurement).
Test: emitted rows carry event_component_status labels; no silent key
shadowing (duplicate-key lint via ast walk in the test).

**B4-11 - CONFIRMED products vanished under junction QC.**
Problem: external CONFIRMED fusion products blocked by junction QC were
dropped without a trace - the one class of evidence an operator most needs
to see when it is excluded.
Fix: excluded products get a ledger (sv_peptide_exclusions.tsv) with reasons
and a provenance count of how many were blocked.
Test: blocked CONFIRMED product appears in the ledger with its reason.

**B4-12 - SV scorer rejected the main branch's auxiliary predictors.**
Problem: run-sv-score accepted only core predictors, so NetMHCstabpan/NetChop
results available in the branch could not be scored together with SV
candidates (the review's closure-evidence item 5: promised aux role, no fake
scores).
Fix: sv-score / sv-run-full accept --netmhcstabpan / --netchop, normalized
into the shared evidence schema; the provenance payload records
predictor_evidence_status so aux completeness is visible.
Test: aux flags parse into schema columns + status recorded; absence does
not fabricate scores.

**B4-13 - adjacent same-strand genes counted as junction support.**
Problem: projected partner-junction matching accepted matches to a same-strand
neighbouring gene - a read spanning a housekeeping neighbour's junction could
support a fusion peptide for the gene next door.
Fix: such matches are relabeled ..._ADJACENT_GENES_REVIEW, their support
reads are zeroed for evidence purposes, and peptides from that support are
blocked; the signal survives as a review flag.
Test: adjacent-gene junction yields the review label, zero counted support,
no peptide.

**B4-15 - WES capture tier HIGH on one-sided evidence.**
Problem: the joint capture interpretability tier granted HIGH when one side
of the event ranked well, letting half the evidence carry the label.
Fix: HIGH requires both sides rank >= 2; one ON-side plus one OFF-side caps
at LOW.
Test: matrix over both-side/one-side/no-side rank combinations.

**B4-17 - capture annotation attached by setattr.**
Problem: SVCluster.capture_annotation was bolted on with setattr - invisible
to the dataclass schema, so any reconstruction of the object dropped it.
Fix: declared typed SVCluster field.
Test: dataclasses.fields() contains the annotation field; round-trip keeps it.

**B4-18 - gene-pair junction keys boosted read support.**
Problem: the SV evidence layer looked up junction reads by gene-pair key when
the exact key missed, crediting a fusion peptide with reads from a different
breakpoint in the same gene pair (the evidence-layer twin of C6-03).
Fix: boosts flow through exact event/coordinate keys only; gene-pair reads
remain in an exploratory column, never in support.
Test: gene-pair-only evidence gives exploratory count but zero support boost.

---

### Package G-e - ranking algebra, endpoint separation, identity vocabularies (commit afa46dd, 2026-11-06 19:53)

Files: src/neoag/evidence_states.py, scoring.py, evidence_consensus.py,
candidate_identity.py, vaccine_events.py, status_words.py, profiles/*.toml;
tests: tests/test_gpackage_ge_20261106.py (33 tests); one pin updated in
tests/test_score_appm_modifiers (strictest-wins cap merge).

**D1-01 - a binding (BA) rank scored in the expression-like (EL) slot.**
Problem: derive_presentation_consensus fell back to binding_rank when EL
fields were missing, so "binds well" was counted as "presented well" - two
different endpoints of different predictors, exactly the mixup the review's
D1 section warns about.
Fix: BA-only evidence stays PRESENTATION_UNASSESSED with reason code
PRESENTATION_EL_UNASSESSED_BA_OBSERVED; BA is recorded as noted, never
scored against EL thresholds.
Test: BA-only row yields UNASSESSED + the reason code (33-test file).

**D1-02 - the same split proven end to end.**
Problem: the review asked not just for the state fix but proof that measured
EL numbers still discriminate.
Fix/test: end-to-end rows through score_evidence_consensus: missing
presentation evidence stays UNASSESSED while a real 99.0/0.0 measurement
scores PRESENTATION_WEAK / R4. Both directions pinned.

**D2-01 - a dead merge function lifted stricter caps when called.**
Problem: merge_peptide_escape overwrote priority_cap with the incoming value
instead of keeping the stricter one - dead code today, a cap-lifting gun
loaded and left on the table.
Fix: deleted; the live merge_escape_flags is strictest-wins and now pinned
by tests (the appm-modifier test pin was updated to the strict semantics).
Test: order pairs through merge_escape_flags keep the stricter cap; the dead
function no longer exists.

**D2-02 - several ad-hoc cap orderings, no total order.**
Problem: cap comparisons were spread across modules as ad-hoc string/letter
checks, and C vs C_CAUTION were treated inconsistently, so two caps applied
in sequence could depend on application order.
Fix: one typed ladder D < C_CAUTION < C < B_CAUTION < B < A shared by all
cap operations; sequential caps are order-independent; off-ladder cap values
are ignored explicitly rather than compared as strings.
Test: pairwise lattice checks + order-permutation of sequential caps.

**D2-03 - ranking cutoffs hardcoded in the function.**
Problem: priority() thresholds (0.75/0.55/0.35) were code constants -
un-versioned ranking policy that can silently change with a refactor.
Fix: thresholds move to [priority_thresholds] in profiles, version-stamped;
the built-in defaults reproduce 0.75/0.55/0.35 exactly (no behaviour change);
every peptide row carries priority_policy_version so results state which
policy made them.
Test: default profile reproduces old boundaries; override changes
assignments; a boundary-ablation test demonstrates sensitivity at each
cutoff.

**D2-04 - orphan peptide rows disappeared from scoring silently.**
Problem: peptides in the ranking input that matched nothing in the raw table
were dropped with no record.
Fix: score() always writes ranked_peptides.unmatched_peptides.tsv - a
header-only file is the visible statement "zero orphans" - classifying each
orphan as absent-from-raw vs filtered-by-enabled_sources.
Test: orphan row lands in the sidecar with its classification; clean input
yields the header-only file.

**D3-02 - the EL tiebreak column mixed in BA rank.**
Problem: the consensus rank tiebreak read a column that could hold a BA rank
when EL was missing (D1-01's ordering twin), and treated missing like zero.
Fix: the tiebreak reads EL fields only; missing values sort last, behind
measured zeros.
Test: BA-rank value cannot act as EL tiebreak; missing sorts behind 0.0.

**D3-03 - contradictory normal-junction evidence graded "safe".**
Problem: the pareto grade summed three normal-junction sources, so
EXACT_MATCH (seen in normal) plus NOT_DETECTED from two others still reached
grade 3 - contradiction averaged into safety.
Fix: the three sources are evaluated separately; the EXACT_MATCH +
NOT_DETECTED contradiction is now grade 1.
Test: the contradiction fixture grades 1, full agreement stays 3.

**D3-04 - the legacy comparator trusted row order.**
Problem: the weighted-rank comparison scored Spearman/Top-N agreement against
whatever row order the baseline file happened to have - a shuffled baseline
would show agreement or disagreement that meant nothing.
Fix: input row order must pass _legacy_baseline_order_status (monotonic
declared rank column) before Spearman/Top-N are emitted; otherwise
NOT_COMPARABLE rows plus a validation row. Baseline source, digest, MISSING
state and order status all land in the run manifest.
Test: unordered baseline -> NOT_COMPARABLE, no correlation numbers; ordered
baseline validates and compares; manifest records provenance.

**D3-05 - event group keys ignored the sample.**
Problem: in a multi-sample input the event/phase group key could be identical
for two samples' same-coordinate events, merging them into one event.
Fix: sample_id is part of the group key when present (the legacy key is
preserved for single-sample rows so single-sample output is unchanged); the
summary gained distinct_sample_ids.
Test: two samples at one coordinate stay two events; single-sample key
byte-identical to before.

**D4-02 - vaccine grouping split the same allele across notations.**
Problem: vaccine peptide-HLA grouping keyed on raw HLA text, so A*02:01 and
HLA-A*02:01 entries for one allele were two groups, each half-populated.
Fix: grouping uses canonical normalize_hla; detail rows keep the original
representation plus a new hla_allele_canonical column so nothing is lost.
Test: mixed-notation rows form one group; originals preserved alongside the
canonical column.

**D4-04 - identity-less rows collided on one digest.**
Problem: rows with neither an id nor coordinates hashed the empty string, so
unrelated unresolved rows shared one identity and could collapse into each
other.
Fix: such rows get row-specific UNRESOLVED_EVENT_ROW digests instead of a
single shared hash.
Test: two identity-less rows yield two distinct digests.

**D4-05 - several incompatible "this cell is empty" vocabularies.**
Problem: candidate_identity, evidence_consensus and vaccine_events each kept
their own list of null texts ("", "NA", "None", "UNASSESSED"...), so a value
empty to one module was a value to another.
Fix: one shared vocabulary - status_words.first_text / is_null_text - used by
all three.
Test: cross-module cases (a "NA" cell is blank to all three paths).

---

### Package G-f - reporting, review and benchmarking honesty (commit 1110ac2, 2026-11-06 20:02)

Files: src/neoag/agent_skills/patient_report.py, reports_dual.py,
report_from_final.py, review.py (open_neo), review_integrity.py,
benchmark_system.py, scoring.py, validation_design.py, evidence_consensus.py,
status_words.py; tests: tests/test_gpackage_gf_20261106.py (38 tests); the
release-audit check-count pin in tests/test_reports_dual.py moved 11 -> 12.
Suite after: 1208 / 119 / 3. Integration subsets: 254 passed.

**E-01 - patient report HTML was an injection surface.**
Problem: agent_skills/patient_report.py interpolated TSV-derived strings
straight into the HTML page, so any field containing markup (a tool comment,
an allele field crafted upstream) could inject active content into the
patient-facing file.
Fix: the HTML build escapes the page title and every section heading and
body via html.escape; the markdown output stays raw markdown (its renderer
does not execute HTML).
Test: a payload value containing <script> and <img onerror> is inert in the
generated page (38-test file).

**E-02 - the release audit's privacy scan missed the deployment roots.**
Problem: release-audit check 10 scanned for private paths using a root list
that did not include the actual deployment roots, so a leaked
/data2/.../patient path in a released report passed the check that exists to
catch exactly that.
Fix: the scan covers /mnt, /root, /home, /data1, /data2, /opt, /srv plus any
absolute two-level root discovered in the run's own provenance, and the
check text now states honestly what is scanned and that the list is a
backstop, not a proof.
Test: a report body containing a /data2/ path fails check 10; a clean report
passes; roots discovered from provenance are covered.

**E-04 (and D3-01) - an all-R3 run looked fine at summary level.**
Problem: the review-integrity overview had an all-R4 check but missed the
all-R3 collapse: a run where every candidate is capped (e.g. because the
safety reference assets were absent) produced a normal-looking summary - the
review's own example being the 88-R3-event patient report. D3-01's actual
demand - being able to tell, from summary level, that caps came from missing
evidence rather than measured unfavorable biology - was resolved by this same
mechanism.
Fix: new no_top_grade_concentration check in review_integrity.py: when there
are >= 20 peptides and zero R1/R2, the run gets a WARN (never a BLOCK, per
the agreed disposition) whose detail summarises the cap reasons parsed from
the per-row evidence_grade_cap_reasons column - including the fraction of
capped rows attributable to absent resources and the top-5 cap codes. The
per-row reason codes already existed; what was missing was the rollup, which
is exactly what both E-04 and D3-01 asked for, and the run manifest already
carries the asset data the rollup reads.
Test: a 20-row R3 fixture with CAP_SAFETY_PARTIAL,CAP_NORMAL_JUNCTION_UNASSESSED
reasons produces the WARN with the concentration figures in its detail; one
R1 row silences it. Overall verdict machine verified: WARN maps to PARTIAL,
not BLOCK.

**E-05 - four modules kept four private phasing vocabularies.**
Problem: scoring, evidence_consensus, review priority and the validation
planner each string-matched phasing values their own way (e.g. equality to
"PHASING_REQUIRED"), so an event marked UNPHASED or PHASING_UNRESOLVED could
skip the phasing ladder in some layers while caught in others.
Fix: one shared vocabulary - status_words.phasing_needs_resolution - with
self-firing tokens (UNPHASED, LOW_CONFIDENCE) and context tokens (REQUIRED,
UNRESOLVED, requiring "PHAS" in the value to avoid unrelated REQUIRED
fields). Bare "LOW" is deliberately excluded: haplotype.py legitimately
writes phase_confidence="low" for every resolved trans pair
(PHASED_NOT_CIS), and the consensus row joins that field in, so a bare LOW
would have put every resolved trans event onto the phasing ladder - a
regression this design avoids and the tests pin. Wired into the D3 cap,
the review priority gate, the scoring cap and note, and the validation
planner.
Test: parametrized helper matrix (fires on PHASING_REQUIRED, UNPHASED,
PHASE_UNRESOLVED, LOW_CONFIDENCE standalone and joined forms; silent on
PHASED_CIS, PHASED_NOT_CIS + low, PHASING_NOT_REQUIRED, negated forms,
REPAIR_REQUIRED) plus a cross-stage matrix showing the same verdict in all
four consumers.

**E-06 - missing purity QC claimed "ASSESSED".**
Problem: report_from_final._purity_qc_status returned ASSESSED whenever the
code path ran, even when no QC artifact existed - a released report could
state the purity was QC'd when no QC was provided.
Fix: absence of a QC artifact returns QC_NOT_PROVIDED; a recognised tool
marker whose status line is absent returns QC_MARKER_WITHOUT_STATUS; only a
genuine status line yields an assessed verdict.
Test: each scenario's exact status string (no marker, marker without
status, PURPLE and FACETS status lines, foreign tool).

**E-07 - a failed profile load silently ran on defaults, and release audit said nothing.**
Problem: if the requested analysis profile failed to load, the report was
built on built-in defaults with no trace in provenance - the report's numbers
stated a configuration that was not used.
Fix: the load path records profile_status LOADED (or FALLBACK_DEFAULTS plus
the error text, with the stub profile labelled) in provenance; the patient
release audit gained check 12 "analysis profile loaded as requested, no
silent fallback", which fails the release on FALLBACK_DEFAULTS and passes on
LOADED or NOT_RECORDED (honest for older artifacts). Release audit is now
12 checks.
Test: bundle-level check 12 pass/fail/NOT_RECORDED; an end-to-end rebuild
with a forced load failure writes FALLBACK_DEFAULTS into both provenance.json
and the audit json, and the LOADED path records LOADED.

**E-08 - a bundle of small reporting honesty fixes (compound, adopted whole).**
Problem: six separate misrepresentations in reports and benchmarks:
(a) multi-tool purity fallback took purities[0] - an arbitrary tool chosen
silently; (b) the release banner was inserted at a fixed line index, so a
layout change put the draft banner in the wrong place; (c) tumor-normal
pairing inferred from file names was reported as if confirmed; (d) the
review's fusion section was spliced into patient HTML by replacing the exact
string "</body></html>", failing silently on any other markup;
(e) the ligandome join overwrote measured 0.0 scores with blanks/fill, and
the APPM multiplier ablation rewrote missing and non-positive multipliers via
`or 1.0` and clamped with min(1.0, ...), manufacturing ablation deltas;
(f) benchmark scoring filled missing tool scores with 0.0.
Fix: (a) median across tools with status
MULTI_TOOL_ARBITRATION_REQUIRED when purities disagree (range > 0.10),
labelled as review guidance that does not replace human arbitration;
(b) the banner anchors to the first section element (fallback: file top);
(c) pairing rows carry pairing_evidence_level
(NORMAL_DNA_DEPTH_EVIDENCE / NO_PAIRING_EVIDENCE /
FILENAME_INFERENCE_NOT_FINGERPRINT) and the wording says "filename
inference, not fingerprint or BAM-Matcher confirmed"; (d) a pure helper
_place_fusion_section inserts before </body>, else </html>, else appends,
and returns the placement label for the artifact; (e) the join preserves
measured scores including 0.0 and fills only blanks; the APPM ablation keeps
a missing multiplier neutral (uncounted), counts observed non-positive ones
in the interpretation instead of rewriting them, and the 1.0 clamp is
removed so the delta is arithmetic; (f) documented: the benchmark's score
matrix keeps its blank fill but labelled rows carry the real score or "", and
the fill is no longer applied over measured values in the join path.
Test: each sub-item unit-pinned in tests/test_gpackage_gf_20261106.py -
median/arbitration values, banner position on varied markup, pairing labels,
placement labels for four markup variants, join zero-keep/fill, multiplier
interpretation text for one non-positive row.

---

## 2. Scope decisions from the adjudication

The campaign implemented the adjudication document (STAGE_REVIEW_ADJUDICATION_
2026-11-06.md) as written. Where the adjudication modified the review, the
modified form is what was implemented; the two items the adjudication closed
outright are recorded in section 3. C4-02 was fixed in this campaign (package
D) per the approved rebuttal's revised scope - the penalty-inversion fix
above is that fix. No item was silently narrowed: each entry above states
what was adopted and, where the review bundled several claims, which parts
landed.

---

## 3. Findings closed before this campaign

**B3a-01 - fusion window coordinate space (closed pre-campaign).**
The review claimed fusion peptide windows were cut without checking that the
breakpoint coordinate and the sliced sequence share a reference frame. The
mechanism was already repaired before this campaign opened, at dev commit
c98fcff: fusion windows ride the breakpoint's own coordinate space
(fusion-protein coordinate preferred with length validation; the explicit
neo breakpoint is the only fallback; invalid rows yield zero peptides and are
QC-flagged in easyfuse_coordinate_qc.tsv), with splice coordinate names
resolved by GTF containment at both the catalog and live normalization
routes. The review's specific "no k-mer check" claim was false against the
dev tree (the proteome k-mer check is applied; see the adjudication's
disproof). Closure evidence: 10 tests at c98fcff plus the live A/B re-run
(30/34 novel peptides) recorded in IMPROVEMENTS.md round 4.

**C1-10 - ranking weight sensitivity (closed pre-campaign).**
The sensitivity knob and bake-off mechanism were already landed (IMPROVEMENTS
entry family 22). The review's residual clause asked for weight re-tuning,
which is an unlabelled ranking change - forbidden under the standing rule
that weight changes require labels and a decided utility model; the
mechanism to explore it exists and is auditable. Closure evidence: the
adjudication document, section on C1-10.

---

## 4. Intentional behavior changes the operator should know

These are consequences of the fixes, not accidents. None of them changed
patient data or ran any pipeline; they change what future runs emit.

Package A: dead classifiers deleted (evidence_state/EvidenceState/STATES,
_status_from_text) - any downstream code that referenced them would now fail
loudly instead of using known-bad matching; no references existed.

Package C: variants whose annotation lacks WT protein context are now
DROPPED (and counted in variants_dropped_no_protein_context) instead of
receiving a fabricated WT peptide - candidate counts can drop on runs where
the WT-protein plugin did not deliver; rows carry wt_peptide_status.

Package D: unresolved clonality multiplies by 1.0 (was 0.65) and
clonality_indeterminate by a new profile knob (default 0.75) - ranking of
clonality-unresolved candidates moves up slightly by design. Blank
self_similarity is now "" + self_similarity_status=UNASSESSED instead of
"0.0" - consumers that parsed "0.0" as a score must read the status column.

Package E: INCOMPLETE junction-read QC defaults to REVIEW (was FAIL-equivalent;
opt-in reject policy available); O3 ORFs now pass the credible-ORF gate;
favorable presentation requires percentile <= 2.0 - rows that previously
earned favorable tiers on predicted values now cap at P1 with a reason.

Package F: records with FILTER="." are excluded by default but ledgered
(parsed/sv_vcf_filter_exclusions.tsv; SV counts may drop slightly versus old
runs, every drop named); orientation-impossible fusion joins are rejected
before peptide generation (fusion peptide counts can drop where the old code
built biologically impossible chimeras); unassessed-normal specificity now
caps at 0.85 (profile knob normal_unassessed_specificity_cap).

Package G-a: unresolved multi-sample VCF roles now abort instead of guessing;
consensus/cross-site keys embed genome build - cross-build runs that
previously merged will now keep separate keys.

Package G-b/G-c: absence is "" not 99.0/0.0/"NA" in predictor and annotation
columns (downstream consumers must tolerate blanks; status columns carry the
meaning); stored tool outputs and composite evidence tables are reused only
with a matching fingerprint sidecar - first run after upgrade re-runs
previously cached steps by design; PRIME composite component is rank-only.

Package G-d: adjacent-gene projected matches no longer count as RNA junction
support (support reads zeroed, peptides blocked - SV peptide counts can
drop); prose/percent purity strings now raise instead of silently parsing.

Package G-e: priority caps share one ladder where C and C_CAUTION merge
(order-independent, strictest-wins); peptides carry a new
priority_policy_version column; score() always writes
ranked_peptides.unmatched_peptides.tsv; vaccine event/phase group keys
include sample_id when present (single-sample output unchanged); vaccine
grouping details gained hla_allele_canonical.

Package G-f: patient release audit is now 12 checks (check 12 = profile
loaded as requested - a silent profile fallback now BLOCKS release); purity
QC vocabulary gained QC_NOT_PROVIDED / QC_MARKER_WITHOUT_STATUS (a missing
QC file no longer reads as ASSESSED, so release audits on such runs now fail
honestly); multi-tool purity disagreement > 0.10 reports
MULTI_TOOL_ARBITRATION_REQUIRED instead of silently taking the first tool;
pairing_status rows gained pairing_evidence_level; patient HTML fusion
section placement is reported as fusion_section_placement.

---

## 5. Status map - all 136 findings

Aggregate: 93 landed in this campaign (packages A-G-f), 2 landed earlier
with closure evidence (B3a-01, C1-10), 1 required no action (positive
citation B3b-15), 4 partially landed (residual named), 36 deferred with
reason. Total 136.

B1 - variant annotation
B1-01  deferred - VEP --pick policy. Changing to a full-CSQ/MANE-preferring
       policy alters which transcript consequences feed every downstream
       variant-to-peptide set; the review itself treats it as "kept on the
       list until an explicit pick-order policy replaces it". An operator
       policy decision with case-level recall impact, not a code-correctness
       fix; B1-04/05 (landed) already guard the WT-column dependency.
B1-02  landed G-a
B1-03  landed G-a
B1-04  landed C
B1-05  landed C

B2 - peptide extraction
B2-01  landed C
B2-02  landed C
B2-03  landed C
B2-04  landed C
B2-05  landed C
B2-06  landed C
B2-07  landed C
B2-08  deferred - verdict "AGREE - PERFORMANCE": proteome k-mer index rebuilt
       each run. No correctness effect; a digest-keyed cache will land when
       memory/time pressure is measured (review's own acceptance: compare
       initial vs cached runs).
B2-09  deferred - verdict "AGREE - DOCUMENTATION": unused codon table and
       labelled minigene back-translation. Documentation/labeling hygiene;
       minigene_nt is already only a construct hypothesis, not used in any
       evidence decision.

B3a - fusion peptide reconstruction
B3a-01 landed earlier - mechanism repaired pre-campaign at dev c98fcff
       (coordinate-space repair, 10 tests, live A/B); review's "no k-mer
       check" claim disproved. Section 3.
B3a-02 deferred - dedicated WT_NOT_APPLICABLE semantics for chimeric
       junctions (verified absent). The review agrees empty paired WT is not
       itself a defect; the residual work is partner-native plus whole-
       proteome comparison design for the fusion evidence layer. The
       generic absence-honesty machinery (C6-01, D4-05, E-statuses) landed
       meanwhile.
B3a-03 deferred - easyfuse event rows still carry hardcoded clonality/
       specificity/expression priors. The labelling pattern exists (G-d
       event_component_status on the SV route); migrating the easyfuse route
       to prior-versus-measured separation is a data-model change for the
       fusion layer, scheduled there.
B3a-04 deferred - BAM-level exact-junction recount for fusion calls; review
       verdict CONDITIONAL (mandatory only per declared evidence grade).
       CALLER_REPORTED_UNVERIFIED is retained and capped as the review's
       fallback clause allows.
B3a-05 deferred - regex-stripping fusion_protein shifts coordinates; the
       review asks to inspect the source format first (some characters are
       benign delimiters). Bundled into the fusion-layer pass with B3a-02/03.
B3a-06 deferred - minigene_nt empty on fusion rows; review: "not inherently
       incorrect", a feature gap only if validation claims nucleotide
       constructs. Same family as B2-09 labeling.

B3b - alternative splicing
B3b-01 landed E      B3b-02 landed E      B3b-03 landed E
B3b-04 landed E      B3b-05 landed E      B3b-06 landed E
B3b-07 landed E      B3b-09 landed E      B3b-10 landed E
B3b-11 landed E      B3b-12 landed E      B3b-13 landed E
B3b-08 deferred - strand-ambiguous ("." ) normal-panel crosscheck; verdict
       CONDITIONAL - SCIENTIFIC REVIEW: the choice between marking such
       crosschecks UNASSESSED or evaluating both strands "under controlled
       ambiguity" is a scientific policy decision the review deliberately
       left open (both options have opposite failure risks). Flagged to the
       operator; current behaviour never calls normal absence from unknown
       strand.
B3b-14 deferred - bundle of seven hygiene claims the review ordered "split
       into separate tickets" (confidence-ladder calibration labels, dead
       branch, HLA-II substring locus match, QC scan count, strict=False on
       secondary sources, splicemutr rejoin build-awareness, manifest
       software_version). The genome-build sub-family landed on the
       consensus/cross-site routes (G-a/C6-06); the splicemutr rejoin site
       itself and the rest stay open as individual tickets.
B3b-15 no action needed - the review's positive citation ("AGREE WITH GOOD
       DESIGN"): INCOMPLETE-never-imputed-as-PASS is retained; the related
       B3b-07 work added the explicit policy treatment it gestures at.
B3b-16 deferred - dead branch in splice2neo; remove-or-implement depends on
       documented branch intent the review itself asks to recover. Hygiene.

B4 - structural variants
B4-01 landed F       B4-02 landed G-d     B4-03 landed G-d
B4-04 deferred - class-II track; verdict "AGREE AS FEATURE GAP" and cross-
       cutting disposition F7: a separately validated future research track
       (allele typing, lengths, predictors, independent safety validation),
       not a correction of class-I output. Current declared scope is
       class-I.
B4-05 landed G-d     B4-06 landed F       B4-07 landed F
B4-08 landed F       B4-09 landed F       B4-10 landed F
B4-11 landed G-d     B4-12 landed G-d     B4-13 landed G-d
B4-14 landed F       B4-15 landed G-d     B4-16 landed F
B4-17 landed G-d     B4-18 landed G-d

B5 - evidence aggregation / schemas
B5-01 landed A       B5-02 landed A       B5-03 landed G-a
B5-04 deferred - adopted in the adjudication but not implemented this
       campaign, recorded honestly: intra-source duplicate keys still
       resolve first-non-empty without a conflict ledger (cross-source
       conflicts were already handled). Severity [M-]; scheduled next to
       the comprehensive_evidence conflict machinery.
B5-05 deferred - compound schema-hygiene bundle the review ordered treated
       "each separately" (drift checks on write_tsv, duplicated field
       registries, tool_run_status versions, mirrored conflict TSVs, purity
       header sniff). Individual tickets; none has a demonstrated wrong-
       output path today.

C1 - HLA presentation prediction
C1-01 landed B       C1-02 landed G-b     C1-03 landed G-b
C1-04 landed G-b     C1-05 landed G-b     C1-06 landed G-b
C1-07 deferred - BA-only allele selection; verdict CONDITIONAL: "defect only
       if the user believes best-predicted means EL/combined". Endpoint
       declaration is an operator-facing choice; landing a default switch
       would silently re-pick alleles on existing cases.
C1-08 landed G-b     C1-09 landed G-b     C1-11 landed G-b
C1-10 landed earlier - sensitivity knob + bake-off already landed
       (IMPROVEMENTS entry 22 family); residual clause = unlabelled weight
       change, forbidden by standing rule. Section 3.
C1-12 deferred - compound ("split into subtests" per review): ic50 "0"
       truthiness, grade() constants, duplicate last-wins, count-only chunk
       resume, empty registry version_args, stale 4.1-vs-4.2 docs. Distinct
       priorities per sub-item; individual tickets.

C2 - immunogenicity composite
C2-01 landed G-c     C2-02 landed B       C2-03 landed G-c
C2-04 landed G-c     C2-06 landed G-c
C2-05 deferred - single cross-layer availability contract (installed / ran /
       output-valid / row-evaluable / synthetic). The worst failure mode it
       describes (stub counted as real evidence) was closed by row-level
       provenance (B) and fingerprint guards (G-b/G-c); the unified contract
       is a design pass across config, runner and reports.

C3 - RNA allele counting
C3-01 landed F       C3-02 landed F       C3-03 landed F
C3-04 partially landed - producer side landed in F (UNASSESSED_COMPLEX kept,
       ERROR rows carry blank numerics); residual: consumers must be shown
       to honour count_status so a stored zero never becomes RNA_NEGATIVE -
       count_status currently has no consumer wiring (verified by grep).
C3-05 deferred - junction VAF semantics in the D1 branch: review verdict
       "AGREE WITH QUALIFICATION" asks for junction-read thresholds/PSI
       instead of SNV VAF under a documented model - a modelling item for
       splice evidence grading.

C4 - cancer cell fraction
C4-01 landed G-d     C4-02 landed D       C4-03 landed G-d
C4-04 landed G-d
C4-05 deferred - joint uncertainty over multiplicity/CN/purity; review:
       report conditional-on-chosen-m or marginalise with simulation
       validation. A statistical-model task (P3 family in the review's own
       repair sequence); the separate sensitivity fields remain.
C4-06 deferred - normal-contamination term; verdict CONDITIONAL: derive the
       exact likelihood and purity definition first - an arbitrary correction
       could double-count. Scientific derivation required before code.
C4-07 deferred - compound: mixed-input QC row sample_id mislabel, external-
       only cluster design note, "NA" sentinel. Review: split tickets.

C5 - HLA/APPM/immune escape
C5-01 landed G-c     C5-02 landed G-c     C5-05 landed G-c
C5-03 deferred - gene-integrity "not_assessed" scores neutral 1.0 but
       displays as intact; fix is a display-vocabulary change (UNASSESSED
       label + justified neutral multiplier). Small, appm display pass.
C5-04 deferred - gene-set drift and PTPN2/SOCS1/IRF9 coverage; verdict:
       needs direction/evidence annotations per gene and curation review -
       blanket penalties would be wrong for some genes (review's own words).
       Research curation, F6 disposition applies.
C5-06 deferred - damaging-variant-plus-expression conflict rule
       over-downgrades true LoF; the fix direction (allele-specific
       expression where appropriate, stop assuming LoF implies low TPM)
       changes scoring semantics and wants the C4/C5 scientific pass.
C5-07 deferred - dead module immune_escape_resistance: today unregistered
       (zero production imports), so no live risk; delete-or-repair decision
       and the duplicated B2M/IFNG constants reconciliation pending.
       Recorded honestly as not-done.
C5-08 deferred - compound: two overall-status formulas disagreeing, path
       truthiness vs existence, completeness denominators; split tickets per
       review.

C6 - peptide safety
C6-01 landed D       C6-06 landed G-a
C6-02 partially landed - per-row cap reasons existed and the run-level
       rollup with absent-resource fraction landed (G-f/E-04, and the D2-04
       orphan ledger); residual: a pre-run missing-assets preflight warning
       is not implemented.
C6-03 deferred - safety-relevant, flagged: the gene-pair normal-junction
       fallback still feeds SAFETY_REJECT (re-verified at
       src/neoag/peptide_safety_gate.py L493-509 this session). Downgrading
       gene-pair hits REJECT->REVIEW loosens a fail-closed safety gate on
       real cases; per the standing rule that safety-behaviour relaxations
       are operator decisions (and consistent with the review's
       "preserve fail-closed policy where intended"), this changes what gets
       rejected and is held for an explicit operator call rather than being
       bundled silently into a hygiene commit.
C6-04 deferred - typed missingness in _event_normal_support (n_alt>=3 with
       missing VAF accepted via default-0.0); same safety-gate behaviour
       family as C6-03, held with it.
C6-05 deferred - hardcoded P2/anchor residues across alleles/lengths; verdict
       "AGREE AS HEURISTIC LIMIT" - allele/length-specific anchor modelling
       is research; the label-only alternative awaits the same gate policy
       decision as C6-03.
C6-07 partially landed - the genome-build sub-claim (cross-site keys from
       the manifest) landed via G-a/C6-06; residual sub-claims (mixed status
       families, hardcoded safety multipliers, HSPC tpm>0 positive evidence,
       underscore gene-pair split) are individual hygiene tickets per the
       review's own "separate tests for each" instruction.

D1 - evidence-state inference
D1-01 landed G-e     D1-02 landed G-e     D1-03 landed A
D1-04 landed A       D1-05 landed A
D1-06 deferred - minimal vs strong junction reads both grade 3; review asks
       for a strength model (depth, uniqueness, mapping quality, caller
       independence, long-read separate) - grading design, not a token fix.
D1-07 deferred - compound: missing rna_vaf as zero (modelled together with
       C3-05), unframed-indel track policy (review: a documented policy
       choice), l3_clonality_score-as-CCF provenance. Split per review.

D2 - weighted ranking
D2-01 landed G-e     D2-02 landed G-e     D2-03 landed G-e
D2-04 landed G-e
D2-05 deferred - composite: exact-0.0 composite triggering the legacy
       fallback, and the sample-level APPM fallback keyed on the allele
       string. The missing-versus-zero machinery exists now (status_words);
       applying it at these two boundary-value sites changes ranking at
       score==0 and is grouped with the ranking-missingness pass. The
       benchmark-side zero handling (same family) did land with E-08.

D3 - evidence consensus
D3-02 landed G-e     D3-03 landed G-e     D3-04 landed G-e
D3-05 landed G-e
D3-01 landed G-f (via E-04) - summary-level cap-reason rollup plus the
       pre-existing per-row evidence_grade_cap_reasons and run-manifest data
       give exactly the missing-evidence-versus-measured-unfavorable
       distinction the finding demanded at summary level.
D3-06 deferred - compound (review: "split defects into tickets"): fragile-
       but-correct grade/cap round-trip, unassessed-counted-zero in the
       exposed numeric (completeness axis partially mitigates), SAFETY_PASS
       with no completeness text assumed complete, event_deduplicate=false
       writing an empty events file, substring manual-review matching
       (over-matches toward MORE review, not less).

D4 - event-level selection
D4-02 landed G-e     D4-04 landed G-e     D4-05 landed G-e
D4-01 deferred - two-stage vaccine artifact contract; the review itself
       states the design is intentional ("documented here") and asks for
       explicit artifact staging/versioned filenames/reviewed-output
       pointers - a workflow-contract/documentation item, not a code defect.
D4-03 partially landed - canonical grouping landed with preservation of the
       original representation and hla_allele_canonical (D4-02); residual:
       a separately declared four-field typed identity versus the effective
       binding-equivalence group (the 3-field collapse itself remains the
       declared equivalence policy).

E - reporting and release
E-01 landed G-f      E-02 landed G-f      E-04 landed G-f
E-05 landed G-f      E-06 landed G-f      E-07 landed G-f
E-08 landed G-f
E-03 deferred - event-level hard-fail BLOCK scoping: distinguishing
       event-global codes (e.g. identity/authenticity failures) from
       peptide-specific ones is a review-level policy classification the
       review itself asks to define before changing ("some hard failures,
       however, are event-wide and should propagate"). Held for that
       classification; today the check stays conservative (blocks).
E-09 deferred - compound ("track each subissue separately"):
       candidate_top_n names the event-row cap, mixed R1-R4/legacy A-D
       labels in one column, APPM section rebuilt by parsing generated HTML,
       benchmark cap-echo assertion, VCF path scraped from the command
       string. Individual tickets; each is small, together they need a
       reporting-input pass.

---

## 6. Verification record

Full pytest after every work package (command:
PYTHONPATH=src python -m pytest from the dev tree). Where the package commit
message recorded the totals: A 962, B 968, C 980, F 1048, G-f 1208 passed;
intermediate packages each ran the full suite green apart from these three.
Known failures throughout: exactly the same three tests, all pre-existing and
untouched by the campaign (verified: identical tracebacks at the pre-campaign
baseline 0dc5391^ and at HEAD). Two are environment-dependent:
test_easyquant_wrapper runs scripts/run_easyquant_sample.sh against bp_quant
and dies with exit 3 on this host's wrapper/binary state;
test_tools_env_config fails in the working tree only because the uncommitted
deployment edit to conf/tools.env.sh redirects NEOAG_PROJECT_ROOT to the
deployed tree (it passes in a clean worktree of the same commit). The third,
test_candidate_identity.py::
test_four_level_identity_keeps_events_but_collapses_same_peptide_hla, is
deterministic, not environmental: the August fixture (05e7b34) expects the
patient-report representative to collapse two fusion rows by normalized
peptide-HLA identity, but a later pre-campaign redesign of
_patient_display_candidate_key (7474afc, 4a22add) keys Fusion display rows by
the :: gene pair precisely to avoid collapsing technical events, and the
minimal fixture has no gene field so both rows fall through to distinct
event-ID keys. Production fusion rows carry gene pairs and do collapse.
The campaign diff never touches either function.

Baseline before the campaign: 954 passed / 119 skipped / 3 known.
Final after G-f: 1208 passed / 119 skipped / 3 known.

Dedicated test files added by the campaign (counts collected from the tree
at report-writing time): test_status_words_20261106 (9),
test_row_provenance_20261106 (6), test_extract_bpackage_20261106 (11),
test_ccf_selfsim_dpackage_20261106 (13), test_splice_chain_epackage_20261106
(27), test_sv_rna_fpackage_20261106 (22), test_rna_allele_counts_script
(rewritten, 8), test_gpackage_ga_20261106 (12), test_gpackage_gb_20261106
(16), test_gpackage_gc_20261106 (24), test_gpackage_gd_20261106 (35),
test_gpackage_ge_20261106 (33), test_gpackage_gf_20261106 (38) - 254 tests
across thirteen files, plus updated pins in pre-existing files
(test_reports_dual check count 11->12, test_score_appm_modifiers cap merge,
test_sv_score stub-binding fix, a few fixture updates).

Integration subset run after G-f with --run-integration (audit_followup,
source_chain, evidence_consensus, open_neo_macro_skills,
rna_evidence_scoring, fusion_consensus, reports_parallel_consensus and the
G-f file): 254 passed.

The test_sv_score.py binding failure fixed in G-d was proven pre-existing at
the previous HEAD by running it in a clean worktree before the fix, so it is
recorded as a side-fix, not a regression introduced here.

Control tree (/data2/caojingqi/neoag_control) was not touched by any package;
all twelve commits are on the dev tree's main branch. The fourteen operator
deployment edits were never staged; each package commit lists only its own
files.

## 7. Round-2 follow-up campaign (2026-11-07)

The follow-up document NEOAG_FIX_CAMPAIGN_OPEN_CONCERNS_HANDOFF_2026-11-07.md
re-examined the campaign's own remedies (ten scrutiny items, a re-examination
table, and a deferral list). The rebuttals for the disputed items were issued
the same day in HANDOFF_RESPONSE_2026-11-07.md; everything the campaign did
not rebut was then implemented in the dev tree, plus two specifically
requested end-to-end tests (B4-10 and B3b-02). Nine commits landed in the
order below (971bcb8 .. c00d220); the full suite after each, and the final
state, are in the verification note at the end of this section.

---

### C3-04 - RNA count statuses now survive the adapter (commit 971bcb8)

Files: src/neoag/adapters/rna_vaf.py, src/neoag/evidence_layer.py,
src/neoag/evidence_states.py; tests: tests/test_evidence_consensus.py
(+1 test).

Problem: the allele-counting script emits per-variant count_status values
(ERROR, CONTIG_NOT_IN_BAM, UNASSESSED_COMPLEX) that describe a failed or
skipped assessment, but the rna_vaf adapter dropped that field, so the
consumer saw counts indistinguishable from a measured site: a pileup that
crashed on one variant was read as RNA_NEGATIVE - "the site was read and no
mutant RNA was seen" - or the legacy status string was trusted over the
numbers. In plain terms, "we could not look" was laundered into "we looked
and saw nothing".
Fix: RNAVafSupport carries count_status/count_error through as
rna_count_status / rna_count_error; evidence_layer's _rna_support_status and
evidence_states' derive_rna_support treat the three technical states as
unassessed in both directions (they also block the status-string-only
confirmation path), while junction support - which does not depend on the
failed ALT pileup - still counts.
Test: test_c3_04_technical_rna_count_states_never_read_as_measured pins all
four readings (technical negative, blocked string confirmation, blocked
numbers, clean counts still confirmed) plus the two _rna_support_status
boundary rules.

---

### B5-04 + D2-02 + D2-05 - join disclosure and two fallback corrections (commit abba2cf)

Files: src/neoag/comprehensive_evidence.py, src/neoag/scoring.py; tests:
tests/test_gpackage_ge_20261106.py (+1 test), tests/test_score_appm_modifiers.py
(+2 tests).

B5-04. Problem: when one source file carried two different values for one
key (two runs, two ORFs, two allele spellings), the join index kept the
first non-empty value and dropped the rest silently. Fix: the merge choice
is unchanged (first value wins), but every kept/dropped pair is now written
to a conflict ledger (INTRA_DUPLICATE_FIELDS columns: source, key_field,
key_value, field, kept_value, dropped_value), so a hidden ambiguity becomes
a visible row.
D2-02. Problem: a priority-cap token outside the A..D ladder was ignored
rather than applied - the permissive reading of a likely schema mismatch.
Fix: an off-ladder cap token now holds the candidate at D (the most
restrictive rung) and names itself in priority_cap_status as
OFF_LADDER_CAP:<LABEL> via unknown_priority_cap_label(); a renamed producer
vocabulary now fails loudly instead of silently un-capping.
D2-05. Problem: the immunogenicity fallback keyed on magnitude, so a
resolved composite of exactly 0.0 was treated as "no data" and a legacy
field could resurrect it. Fix: availability decides - the fallback runs only
when immunogenicity_source is "none"; a measured 0.0 stays authoritative.
The sample-summary HLA-LOH penalty additionally names its provenance
(appm_multiplier_reason starts sample_level_appm_summary+hla_loh_allele_string:),
distinguishing it from a per-candidate APPM modifier row.
Test: one cap-algebra test and two scoring tests in the named files.

---

### C3-01 - genome-build fingerprint preflight (commit 180385c)

Files: scripts/rna_allele_counts_pysam.py; tests:
tests/test_rna_allele_counts_script.py (8 -> 11 tests).

Problem: the BAM-vs-reference preflight checked contig names only. GRCh37
and GRCh38 use the same chr1..22/X/Y spellings, so a 38-labelled BAM over a
37 reference (or vice versa) passed preflight and produced silently wrong
coordinates. Fix: GENOME_BUILD_FINGERPRINTS holds exact primary-assembly
contig lengths for GRCh37 and GRCh38; the preflight votes on observed BAM
lengths and classifies build_check= as BUILD_CONFIRMED / BUILD_MISMATCH /
BUILD_INFERRED / BUILD_UNASSESSED; --expected-build (aliases hg19/b37/grch37,
hg38/b38/grch38) is enforced and a confirmed mismatch aborts unless the
audited --allow-build-mismatch override is given; the status also lands in
the final wrote-line. An ambiguous or novel assembly (e.g. T2T-CHM13 chr1
length 207734681, which matches neither fingerprint) is reported as no-match
rather than guessed.
Test: fingerprint voting (exact, ambiguous mix, unknown assembly), honest
inference paths, and the abort-on-confirmed-mismatch path (3 new tests).

---

### C6-03 - junction match levels: exact coordinates vs gene-pair context (commit ac7f1d3)

Files: src/neoag/peptide_safety_gate.py, src/neoag/schemas.py; tests:
tests/test_safety_completeness.py (+2 tests).

Problem: the normal-junction lookup keys candidates by their exact junction
coordinates and - as a fallback - by the event's gene pair. A hit through the
gene-pair fallback alone (which only says the gene pair is junction-capable
in normal tissue) drove the same hard reject as an exact coordinate match,
letting context stand in for proof. Fix: both the peptide row and the event
rollup now carry normal_junction_match_level = EXACT_NORMAL_JUNCTION_MATCH /
GENE_PAIR_ONLY_CONTEXT; the reject reason distinguishes the two
(normal_junction_seen[_gene_pair_only], low_level_ variants); the profile
knob safety.normal_junction_gene_pair_only_action = reject (default - keeps
the historical production behavior, now honestly labelled) or review (caps
the context-only case at CAUTION/review). The new columns are declared in
PEPTIDE_SAFETY_FIELDS / EVENT_SAFETY_FIELDS here and in schemas.py.
Test: gene-pair-only case under both knob settings, and the exact-coordinate
case keeping its own level and reason under both.

---

### B4-10 + B4-07 - SV cluster evidence restricted to representative-compatible members (commit 65a7521) - requested test 1

Files: src/neoag/sv/sv_merge.py, src/neoag/sv/sv_filter.py,
src/neoag/sv/phase1.py, src/neoag/sv/schemas_sv.py; tests:
tests/test_sv_rna_fpackage_20261106.py (22 -> 23 tests).

Problem: SV clustering is a transitive closure, so a chain
A~B~C can put C in the cluster even when C is incompatible with the
representative A (opposite strand, or breakpoints beyond the cluster
distance). Every aggregate then mixed incompatible evidence: VAF, caller
agreement, split-read/precision support, the sniffles channel and the normal
penalty. Decision requested: should only representative-compatible evidence
count? Yes - implemented that way, symmetric across all channels.
Fix: SVCluster keeps the permutation-invariant union-find clustering, but
aggregation uses only records pairwise _compatible with the representative
under the cluster's own rule (same type family, both canonical endpoints
within the distance, equal strand where both known). New disclosure columns
keep the exclusions visible: cluster_record_count, compatible_record_count,
cluster_all_callers, excluded_incompatible_record_ids,
excluded_incompatible_tumor_alt_max, and event_evidence_basis =
REPRESENTATIVE_COMPATIBLE@v1 (exclusions present) vs ALL_MEMBERS_COMPATIBLE;
the confidence scorer gains the reason token
incompatible_cluster_members_excluded.
Demonstration (the requested chain, one permutation-invariant cluster):
A = GRIDSS +/+ 100000-200000 alt 40/depth 50, B = Manta ././ 100150-200120
alt 2/depth 30 (compatible), C = SvABA -/- 100300-200300 alt 200/depth 210
(incompatible strand and distance). Representative A. VAF is now
(40+2)/(50+30) = 0.8000 instead of (40+2+200)/(50+30+210) = 0.952381; the
event-confidence score falls from 0.85 to 0.75 (caller agreement 2/2
compatible instead of C's third caller counted); the downstream clonality
and tumor-specificity rollup reads the corrected row (1.0 / 0.8500). C stays
disclosed (excluded record id, its alt max 200, basis label).
B4-07: the tumor-specificity rules are emitted as a versioned policy string
specificity_rules@v1:... per raw SV event so the scoring contract is in the
data, not only in the code.

---

### B3b-02 - merged splice statistics stay a coherent donor tuple (commit b69657e) - requested test 2

Files: src/neoag/splice/pipeline.py, src/neoag/splice/schemas.py,
src/neoag/splice/adapters/irfinder.py, src/neoag/splice/consensus.py; tests:
tests/test_splice_chain_epackage_20261106.py (27 -> 31 tests),
tests/test_splice_v050_layer.py (vocabulary pins updated).

Problem: rows of one event were merged field-by-field with max-absolute
values, so {psi, delta_psi, qvalue} could come from three different rows - a
statistical statement no caller ever produced. Requested test: two rows with
different delta-PSI/q-value pairs must merge to a biologically coherent,
exactly-sourced result.
Fix: the statistics fields {psi, delta_psi, qvalue, outlier_score} are taken
whole from one donor row, chosen by min over the rank tuple
(qvalue or +inf, -|delta_psi|, -|psi|, source_record_ids) - rule
MIN_QVALUE_THEN_MAXABS_DELTA_PSI@v1 (smallest q-value wins, strongest effect
breaks ties, deterministic id tiebreak). The merged row names the donor
(stats_source_record_ids) and the rule (stats_selection_rule); when the old
per-field merge would have produced a different (cross-paired) result, a
conflicts-table entry STATS_TUPLE_CROSS_PAIRING_AVOIDED (severity WARNING,
action COHERENT_TUPLE_SELECTED) documents what was avoided. Identical rows
merge without any conflict.
Demonstration (the requested pair): rows chr1:151-200:+ GENE1 with
(psi 0.5, delta_psi -0.9, q 0.4) and (psi 0.6, delta_psi 0.1, q 0.00001)
merge to the second row's tuple verbatim (q 1e-5 wins) with
stats_source_record_ids equal to that row's source_record_ids; the cross-
pairing avoided entry names the tuple the per-field rule would have built.
Related in the same file: IRFinder's IRratio is vocabulary-separated into
its own ir_ratio column (psi left empty) - an intron-retention ratio is not
a PSI; the B3b-06 guard test pins that INCOMPLETE junction QC never reads as
verified (decision stays REVIEW + EXPLORATION_EVIDENCE_INCOMPLETE under both
policies); B3b-12's FAVORABLE_PRESENTATION_PERCENTILE = 2.0 gained its full
provenance comment (pVACbind EL presentation percentiles per allele x
epitope-length background, companion EL score field, pre-declared
exploratory threshold).

---

### D4-04 + B2-01 - identity derivation labels; loud WT-context loss (commit fc11fcd)

Files: src/neoag/candidate_identity.py, src/neoag/vep/extract_peptides.py;
tests: tests/test_candidate_identity.py (+1 test, fixture corrected),
tests/test_extract_bpackage_20261106.py (+1 test).

D4-04. Every candidate identity row now names how its identity was derived:
event_identity_status = EXACT_PRIMARY_ID / COORDINATE_DERIVED_ID /
UNRESOLVED_ROW_FALLBACK, so a fallback identifier can no longer be mistaken
for a primary one. The third pre-existing suite failure was closed here: the
stale August identity fixture predates the :: display-key redesign (7474afc,
4a22add) and lacked the gene pair real fusion rows carry; with the gene field
added it pins the intended collapse (two rows, one display key) and the
identity status remains green. The known-failure count is therefore 2 from
this commit on.
B2-01. Variants dropped for missing protein context were only a counter.
The extractor now emits a stderr WARNING and a sensitivity_warning summary
key, and the message names the scope: row-level ("see
variants_no_peptides_by_consequence for the per-reason funnel") when the CSQ
has a WildtypeProtein column, versus batch-wide ("the protein VEP plugin
likely did not run ... re-annotate before trusting sensitivity") when the
column is absent - the difference between one broken row and a silently
missing annotation stage. B2-01's never-fabricate-WT behavior itself is
unchanged and stays pinned by its package-B tests.

---

### B4-02 - FASTA reader follows the real htslib .fai convention (commit c8550f9)

Files: src/neoag/sv/reference.py; tests: tests/test_gpackage_gd_20261106.py
(35 -> 42 tests).

Problem: the byte-offset index treated .fai column 4 as the line width
including the newline (chars_per_line = col4 - 1). samtools faidx puts
bases-per-line in column 4 and the byte width in column 5, so against any
real samtools-indexed reference the reader strode one base early from line 2
onwards - every micro-homology and protein-reconstruction fetch would read
shifted sequence, silently. The shipped tests passed because their helper
wrote the same nonstandard convention (col4 = col5 = width + 1); the
HANDOFF's real-reference-compatibility concern was correct.
Fix: the parser now implements the spec (col4 bases per line, col5 = col4 +
1 or + 2 separator width), refuses anything else, and verifies the sidecar
against the file - an over-declared length (file truncated after indexing)
or an under-declared one (bases continuing past the declared end) falls back
to the honest streaming scan instead of fetching garbage. Content beginning
with the gzip magic is never byte-addressed regardless of file name (whole-
file fallback, which decodes correctly); a genuine single-line plain FASTA
(-1 -1 columns) addresses directly; the scanner rejects ragged line layouts
(non-uniform widths or blank lines inside sequence) for byte arithmetic and
falls back honestly. The old 71/71 sidecar convention is now refused rather
than trusted.
Tests: samtools LF and CRLF indices (asserting the parsed stride),
single-line -1/-1, gzip-magic-under-plain-name, stale/malformed sidecar
fallback in both directions, first/last/base-boundary/empty/out-of-range
fetches under indexed mode, ragged-layout fallback, and the no-materialisation
memory guard.

---

### E-03 + E-04 - typed hard-fail scopes; risk-gated no-R1/R2 summary (commit c00d220)

Files: src/neoag/hard_fail_scopes.py (new), src/neoag/evidence_consensus.py,
src/neoag/open_neo/review.py, src/neoag/open_neo/review_integrity.py; tests:
tests/test_open_neo_macro_skills.py (+4 tests).

E-03. Problem: a hard failure on one peptide - exact normal-proteome match,
one lost restricting allele - blocked the whole event, including defect-free
siblings, through the event rollup (any hard-failed peptide made a promoted
event an integrity violation) and the review's representative check.
Conversely, genuinely global failures (sample identity, invalid coordinates)
must keep blocking everything.
Fix: every hard-fail code is typed with one of five scopes - SAMPLE_GLOBAL,
EVENT_GLOBAL, ORF_SPECIFIC, PEPTIDE_SPECIFIC, PEPTIDE_HLA_SPEC; unknown codes
default to EVENT_GLOBAL, so a new reason never loosens safety by omission
(the registry also proves no binder/prediction code is hard-fail at all:
weak presentation caps the grade, it does not invalidate siblings).
Propagation: peptide rows carry hard_failure_scopes (CODE:SCOPE pairs); the
event rollup reports hard_failure_event_blocking, hard_failure_event_scope_
codes, hard_failure_peptide_local_codes and hard_failure_event_block_basis
(EVENT_OR_SAMPLE_SCOPE_CODE / ALL_PEPTIDES_HARD_FAIL / PEPTIDE_LOCAL_ONLY /
NONE - local defects still block when no clean candidate remains). Review
excludes the event only when the rollup (or a rep-level event-scope code)
says so; a defect confined to the representative holds it
(HOLD / REVIEW_REPRESENTATIVE_LOCAL_HARD_FAIL - it may not be promoted, but
the event survives, and select_first_batch only admits it through the single
manual-review lane); a clean event with a defective sibling advances with
the sibling named in the reason (EVENT_SIBLING_LOCAL_HARD_FAIL=...).
review_integrity's hard_fail_not_promoted check follows the typed rollup,
keeps the per-peptide check, counts promoted events with visible sibling
defects (promoted_events_with_sibling_local_hard_fail=N) so non-blocking is
provably disclosed, and keeps the pre-E-03 rollup gate for artifacts without
typed columns. All four required scenarios are tested: exact-normal-match
sibling, invalid event identity (still excludes and integrity-BLOCKS if
promoted anyway), different-HLA-loss, and every-peptide-defective.
E-04. The "no top-grade candidate" summary now warns only when the batch is
actually exposed - no R1/R2 AND (>= 20 peptides, or missing-resource caps in
play, or the top-ranked candidate itself capped) - instead of on every quiet
batch; the cap-completeness counts stay unconditional.

---

### Round-2 verification note

Full suite after every commit (same command as section 6). Final state:
1235 passed / 119 skipped / 2 known failures. The known-failure count dropped
from 3 to 2 because the third failure (the stale identity fixture) was fixed
in fc11fcd with the :: display-key root cause already proven at section 6;
the two remaining failures are the same environment-dependent pair (the
easyquant wrapper against this host's bp_quant state, and the tools.env
external-root test passing in a clean worktree only). Suite growth: 1208 ->
1235 (+27). New/updated pins this round: +1 C3-04 consumer test, +1 cap
test and +2 D2-05 tests, +3 C3-01 tests, +2 C6-03 tests, +1 B4-10 chain
test, +4 splice tests (coherent tuple, identical-row merge, guard, and the
updated v050 ir_ratio pins), +1 identity status test, +1 B2-01 warning test,
+7 B4-02 FASTA boundary tests, +4 E-03 scope tests. Four round-2 changes
(C3-04 consumer branch, D2-05, C6-03, B2-01) initially landed from the
earlier session without pins; those pins were written this round before
landing and are included in the counts above. Targeted integration runs
alongside the full suite: macro-skills 100 passed, consensus + GF + fusion +
parallel-consensus 103 passed, extract + safety set 48 passed (+ the five
appended pins green).

Nine commits: 971bcb8 (C3-04), abba2cf (B5-04/D2-02/D2-05), 180385c (C3-01),
ac7f1d3 (C6-03), 65a7521 (B4-10/B4-07), b69657e (B3b-02), fc11fcd
(D4-04/B2-01), c8550f9 (B4-02), c00d220 (E-03/E-04). The fourteen operator
deployment edits were again never staged; each commit lists only its own
files; the control tree was not touched.

---

This report is mirrored to /data1/neoag/install_logs/REVIEW_FIXES_2026-11-06.md
with an md5 comparison on delivery.

