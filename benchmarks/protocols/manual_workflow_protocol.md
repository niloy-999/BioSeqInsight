# Manual workflow reference protocol

This document defines the "existing tools" baseline that BioSeqInsight is
compared against. It is published for one reason: a comparison is only
meaningful if the thing being compared to is written down. A reader who
thinks this baseline is unfair can point at the specific step they disagree
with.

`benchmarks/scripts/run_workflow_comparison.py` encodes exactly these step
lists. If this document and that script ever disagree, the script is wrong
and should be corrected to match.

---

## Ground rules

* The baseline uses **free, public, web-based tools**, because that is what
  the target users of a desktop teaching-and-research tool actually use. A
  comparison against a fully scripted Biopython pipeline would be a different
  study; see "What this baseline is not" below.
* Each step is one action a user performs: opening a page, pasting text,
  setting a parameter, reading a value, or writing a value down.
* A step is marked as a **manual data transfer** when information moves
  between tools by copy-paste or retyping. These are counted separately
  because they are where transcription errors come from.
* The baseline assumes a competent user who knows which tools to use. Time
  spent searching for a tool is not counted, which favours the baseline.

---

## T1 — Single-sequence characterisation

Goal: GC content, melting temperature, longest ORF, and the molecular weight,
pI and GRAVY of the ORF product.

| # | Step | Application | Manual transfer |
|---|---|---|---|
| 1 | Open a GC-content calculator | web tool A | |
| 2 | Paste the sequence | web tool A | yes |
| 3 | Read and record GC content | notes | yes |
| 4 | Open a melting-temperature calculator | web tool B | |
| 5 | Paste the sequence again | web tool B | yes |
| 6 | Set primer and salt concentrations | web tool B | |
| 7 | Read and record Tm | notes | yes |
| 8 | Open an ORF finder | web tool C | |
| 9 | Paste the sequence again | web tool C | yes |
| 10 | Choose genetic code and minimum ORF length | web tool C | |
| 11 | Read the ORF table | web tool C | |
| 12 | Copy the longest ORF peptide | web tool C | yes |
| 13 | Open a protein parameter tool | web tool D | |
| 14 | Paste the peptide | web tool D | yes |
| 15 | Read molecular weight, pI and GRAVY | web tool D | |
| 16 | Transcribe all values into a notebook or spreadsheet | notes | yes |

**16 steps, 5 applications, 7 manual transfers.**

Note that steps 6 and 10 are parameter choices that each tool defaults
differently, and that nothing in this workflow records what was chosen.

BioSeqInsight equivalent: paste, click *Full analysis*, click *Export CSV*.
**3 steps, 1 application, 1 manual transfer.**

---

## T2 — Structure retrieval with an identity check

Goal: obtain a structure and establish whether it is the queried protein.

| # | Step | Application | Manual transfer |
|---|---|---|---|
| 1 | Open the UniProt entry | UniProt | |
| 2 | Copy the canonical sequence | UniProt | yes |
| 3 | Open the RCSB search page | RCSB PDB | |
| 4 | Paste the sequence and run a sequence search | RCSB PDB | yes |
| 5 | Inspect the hit list and choose an entry | RCSB PDB | |
| 6 | Open the entry page and check its UniProt cross-reference | RCSB PDB | |
| 7 | Download the coordinate file | RCSB PDB | |
| 8 | Extract the chain sequence from the file | text editor | |
| 9 | Open a pairwise alignment tool | EMBOSS / BLAST | |
| 10 | Paste the query sequence | EMBOSS / BLAST | yes |
| 11 | Paste the chain sequence | EMBOSS / BLAST | yes |
| 12 | Run the alignment; read identity and coverage | EMBOSS / BLAST | |
| 13 | Decide whether the match is acceptable | judgement | |
| 14 | If not, open the AlphaFold DB entry | AlphaFold DB | |
| 15 | Download the predicted model | AlphaFold DB | |
| 16 | Read the mean pLDDT | AlphaFold DB | |
| 17 | Record the decision and its evidence | notes | yes |

**17 steps, 7 applications, 5 manual transfers.**

Steps 8–13 are the identity check. They are the steps most often skipped,
and skipping them is precisely the failure that produced v1.0's misleading
"25/25 structures retrieved" result, in which 7 of the 25 were a different
molecule. A baseline that omitted these steps would be a faster baseline and
a dishonest one.

BioSeqInsight equivalent: enter the accession, click *Retrieve and validate*,
read the mapping level. **3 steps, 1 application, 1 manual transfer** — and
the identity check cannot be skipped, because the mapping level is always
shown.

---

## T3 — Batch of fifty sequences

| # | Step | Application | Manual transfer |
|---|---|---|---|
| 1 | Split the FASTA file into individual sequences | text editor | |
| 2 | For each of 50 sequences, repeat T1 (16 steps each) | web tools A–D | yes |
| 3 | Paste each result into a spreadsheet | spreadsheet | yes |
| 4 | Check for transcription errors | spreadsheet | |
| 5 | Format the table | spreadsheet | |

Expanded: **50 × 16 + 4 = 804 steps.**

BioSeqInsight equivalent: **3 steps**, independent of the number of
sequences.

---

## T4 — Reproduce an earlier analysis

| # | Step | Application |
|---|---|---|
| 1 | Find the notebook entry or spreadsheet | notes |
| 2 | Identify which web tools were used | notes |
| 3 | Check whether each still exists and behaves the same | web tools |
| 4 | Guess the parameters that were used | judgement |
| 5 | Repeat the whole analysis | web tools A–D |
| 6 | Compare with the recorded values and reconcile differences | spreadsheet |

**6 steps** — and steps 3 and 4 may be impossible. Whether this task
succeeds at all depends on what the analyst wrote down months earlier.

BioSeqInsight equivalent: open the project, read the stored settings snapshot
and software version, re-run. **3 steps**, and the parameters are recorded
rather than guessed.

---

## T5 — Hand the work to a collaborator

| # | Step | Application |
|---|---|---|
| 1 | Collect the input files | file manager |
| 2 | Collect the downloaded structures | file manager |
| 3 | Collect the results spreadsheet | file manager |
| 4 | Write down which tools and settings were used | notes |
| 5 | Zip the folder by hand | file manager |
| 6 | Send it and explain the layout | email |

**6 steps.** BioSeqInsight equivalent: *Project → Export project as zip*,
then send. **2 steps**, and the archive re-imports into a working project.

---

## What this baseline is *not*

* It is **not** a comparison against a scripted pipeline. A bioinformatician
  who writes a Biopython script will beat both routes for a batch of fifty
  sequences, and should. BioSeqInsight's claim is about users who do not
  write that script — the same audience that a desktop GUI exists for. The
  manuscript must state this scope explicitly rather than implying a general
  superiority claim.
* It is **not** a claim that the individual web tools are inaccurate. They
  are not. The cost being measured is fragmentation: re-entry, unrecorded
  parameters, and an identity check that nothing forces you to perform.
* It is **not** a measurement of human time. Step counts are a proxy.
  Human completion time and error rates require the usability study in
  `usability_study_protocol.md`.

---

## Checking these counts

The step lists live in `run_workflow_comparison.py` as data. To see them,
and the derived counts:

```bash
python benchmarks/scripts/run_workflow_comparison.py
```

To dispute a count, change the step list, re-run, and the totals update. That
is the intended way to disagree with this baseline.
