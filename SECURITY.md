# Security policy

## Supported versions

| Version | Supported |
|---|---|
| 2.0.x | yes |
| 1.0.x | no — archived at the `v1.0-archive` tag for reproducing published results only |

## Reporting a vulnerability

Please report security issues privately, through GitHub's "Report a
vulnerability" form on the Security tab, rather than opening a public issue.

Include what you did, what happened, and what you expected. A minimal
reproduction helps most. You can expect an acknowledgement within a week and
an assessment within two.

## What is in scope

BioSeqInsight is a desktop and command-line client. The areas where a defect
could realistically cause harm are:

* **Archive extraction.** `project import` reads zip archives. Path traversal
  ("zip slip") is explicitly defended against and tested
  (`test_import_refuses_unsafe_paths`). A bypass is a valid report.
* **Generated HTML.** Batch reports and the structure viewer embed
  user-controlled text. Identifiers are HTML-escaped and coordinate data is
  embedded via `json.dumps`. An escaping bypass is a valid report.
* **Cache integrity.** Cached payloads are checksummed; a tampered file is
  discarded rather than served. A way to defeat that check is a valid report.
* **Transport.** The stdlib transport refuses non-HTTPS URLs. A way to induce
  a plaintext request is a valid report.
* **Dependency supply chain.** The package has no required runtime
  dependencies, which is deliberate. A change that silently adds one is worth
  reporting.

## What is out of scope

* The availability or behaviour of RCSB PDB, EMBL-EBI, ESM Atlas or UniProt.
* 3Dmol.js, which is loaded from a CDN in the generated viewer page. Report
  issues in it upstream.
* Scientific disagreement about a calculation. That is a bug report or a
  discussion, not a security issue, and is very welcome as an issue.

## Handling of user data

BioSeqInsight runs locally. It sends a sequence to an external service only
when you ask it to retrieve a structure, and only to the resource being
queried. Logs record a salted hash and the length of a sequence rather than
the sequence itself, so a shared log file does not leak unpublished data.
Offline mode (`--offline`) disables all outbound requests.
