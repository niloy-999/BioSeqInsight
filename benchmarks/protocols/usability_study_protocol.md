# Usability evaluation protocol

**Status: this is a protocol, not results.** No participant data is included
in this repository. The protocol is published here so that the study can be
run, checked and repeated, and so that reviewers can see exactly what will
and will not be claimed from it.

This document specifies the study that supplies the human half of the
"added value" evidence. The machine half — interaction counts and machine
timings — is produced by `benchmarks/scripts/run_workflow_comparison.py` and
needs no participants.

---

## 1. Research questions

**RQ1.** Does an integrated workflow reduce task completion time compared
with the equivalent multi-tool manual workflow?

**RQ2.** Does it reduce transcription and interpretation errors, in
particular the error of accepting a retrieved structure as the queried
protein when it is not?

**RQ3.** How do users rate the usability of the software (SUS), and what
workload do they report (NASA-TLX)?

**RQ4.** Which parts of the interface cause difficulty, and what do users ask
for that is missing?

RQ2 is the one that matters most for this software. A faster route to a wrong
answer is not an improvement.

---

## 2. Design

A **within-subjects, counterbalanced** design. Every participant completes
matched tasks under both conditions:

* **Condition M (manual):** the participant's own choice of web tools,
  following `manual_workflow_protocol.md`.
* **Condition B (BioSeqInsight):** the same tasks in the software.

Within-subjects removes between-participant variance in molecular biology
experience, which would otherwise swamp the effect with any realistic sample
size. Order is counterbalanced (half do M first, half do B first) and the
task sets are matched but not identical, so that the second run cannot be
completed from memory of the first.

**Blinding.** Sessions are facilitated by someone who is not an author of the
software where possible. If that is not possible, the facilitator follows a
written script and must not intervene beyond it; this limitation is reported.

---

## 3. Participants

**Target: 12 participants, minimum 8.**

Recruited from three strata, to check that the result is not an artefact of
one group:

| Stratum | Target n | Definition |
|---|---|---|
| Undergraduate / MSc students in biology or biotechnology | 4 | Has taken a molecular biology course; little or no command-line use |
| PhD students and postdocs in the life sciences | 4 | Uses sequence tools regularly; occasional scripting |
| Bioinformatics practitioners | 4 | Writes analysis code routinely |

Exclusion criterion: prior use of BioSeqInsight.

**On sample size.** Twelve participants will not support a precise effect
size. It is sufficient to detect a large within-subjects difference in
completion time (paired design, α = 0.05, power 0.8, detectable d ≈ 0.9) and
to surface the large majority of usability problems. The report must state
this limitation rather than implying more precision than the design provides.
Report effect sizes with confidence intervals, not only p-values.

---

## 4. Tasks

Matched task pairs, one per condition. Set A and Set B use different but
equivalent inputs.

**Task 1 — Sequence characterisation (target ~10 min manual).**
Given a coding sequence, report GC content, melting temperature under stated
conditions, the longest ORF, and the molecular weight and isoelectric point
of its product.

**Task 2 — Structure retrieval with an identity check (target ~15 min manual).**
Given a UniProt accession, obtain a structure and state, with evidence,
whether it is the queried protein, a homologue, or something else.

*Task 2 includes one deliberate trap:* one accession in each set is a protein
whose top sequence-search hit is a different molecule (the preproinsulin →
insulin-degrading enzyme case from the v1.0 audit). Whether participants
detect this is the primary measure for RQ2.

**Task 3 — Small batch (target ~20 min manual).**
Produce one table of properties for ten sequences supplied as a FASTA file.

**Task 4 — Hand over the work.**
Package everything a colleague would need to verify and continue the
analysis.

Tasks are capped at 25 minutes each. A task that hits the cap is recorded as
censored at 25 minutes and reported as such, not discarded.

---

## 5. Procedure

1. Consent and briefing (5 min).
2. Background questionnaire: field, career stage, tools used weekly, prior
   experience with the resources involved (5 min).
3. Condition 1 tasks, think-aloud, screen and audio recorded (45 min).
4. NASA-TLX for condition 1 (3 min).
5. Break (5 min).
6. Condition 2 tasks (45 min).
7. NASA-TLX for condition 2; SUS for BioSeqInsight (5 min).
8. Semi-structured interview (10 min): what was confusing, what was missing,
   what would stop you using this.

Total: approximately two hours per participant.

**Training.** Condition B begins with a fixed five-minute demonstration of a
task that is *not* in the study set. Condition M begins with a statement of
which tool categories are permitted. Neither condition receives help once the
task starts, beyond re-reading the task description.

---

## 6. Measures

**Primary**

* Task completion time per task per condition (seconds, censored at 1500 s).
* Correctness of the reported values, scored against a pre-computed answer
  key by two independent scorers; disagreements resolved by discussion.
* **Identity errors:** whether the participant accepted a non-matching
  structure as the queried protein. Binary, per task-2 trial.

**Secondary**

* SUS score (0–100) for BioSeqInsight.
* NASA-TLX raw scores for both conditions.
* Number of distinct applications opened (from the screen recording).
* Number of manual copy-paste or retyping operations (from the recording).
* Transcription errors: values that differ from what the tool displayed.

**Qualitative**

* Think-aloud utterances coded for confusion, surprise and satisfaction.
* Interview responses, thematically coded by two coders.

---

## 7. Analysis plan

Fixed before data collection.

* Completion time: paired analysis across conditions. Wilcoxon signed-rank as
  the primary test (no normality assumption at this sample size), with
  medians and bootstrap confidence intervals reported. Paired *t*-test
  reported as a secondary check.
* Identity errors: McNemar's test on the paired binary outcomes.
* SUS: mean with a 95% confidence interval; interpreted against the
  established adjective scale rather than a pass mark.
* NASA-TLX: paired comparison per subscale, with the multiple-comparison
  correction stated in advance (Holm).
* Qualitative: inductive thematic coding; inter-coder agreement reported as
  Cohen's κ.

**Pre-registration.** This plan is to be timestamped before the first session
(e.g. in the repository history or on OSF) and the registration referenced in
the manuscript.

**Stopping rule.** Data collection stops at the target n. No interim
analysis, no adding participants after seeing results.

---

## 8. Ethics and data handling

* Institutional ethics approval must be obtained before the first session.
  The approval reference is to be recorded here and in the manuscript. If the
  institution deems the study exempt, the exemption determination is recorded
  instead.
* Written informed consent, covering screen and audio recording and the reuse
  of anonymised data.
* Participation is voluntary, withdrawal is possible at any point without
  giving a reason, and withdrawn data is destroyed.
* Recordings are stored encrypted, accessible to the named study team only,
  and deleted after the analysis is published.
* Published data is anonymised: participants are identified by stratum and a
  number (for example "P07, postdoc"). No direct quotations that could
  identify an individual or institution.
* The anonymised per-participant measure table is published with the
  manuscript so the analysis can be checked.

---

## 9. Threats to validity, stated in advance

* **Experimenter bias.** The authors built the software. Mitigations: an
  independent facilitator where possible, a written script, pre-registered
  analysis, published raw scores.
* **Novelty effect.** BioSeqInsight is new to participants and the manual
  tools may be familiar, which biases *against* the software for time and
  *towards* it for attention. Reported, not corrected for.
* **Task selection.** The tasks were chosen by the authors and match what the
  software does well. This is why the manual protocol is published: a reader
  can judge whether the tasks are representative.
* **Small sample.** See section 3.
* **Single site.** Participants come from one institution; generalisation is
  limited accordingly.

---

## 10. Reporting

The manuscript must report:

* the achieved sample size and composition;
* medians, interquartile ranges and effect sizes with confidence intervals
  for every primary measure;
* the identity-error rate in both conditions, as the headline usability
  result;
* every task that hit the time cap;
* the SUS score with its confidence interval;
* the limitations in section 9, in the limitations section rather than in a
  footnote.

If the study is not run before submission, the manuscript must say so plainly
and confine its claims to the objective interaction counts and machine
timings, which require no participants.
