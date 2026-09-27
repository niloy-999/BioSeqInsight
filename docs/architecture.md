# Architecture

This document explains where the boundaries in the codebase are and why they
are there. The short version: version 1.0 interleaved computation, I/O and
display in three flat modules, which made the scientific code impossible to
test without a display server. Everything here follows from fixing that.

---

## Layers

```
┌─────────────────────────────────────────────────────┐
│  gui/            cli.py                             │  presentation
├─────────────────────────────────────────────────────┤
│  workflows/      batch runner, projects             │  orchestration
├─────────────────────────────────────────────────────┤
│  structures/     providers, validation, manager     │  domain services
├─────────────────────────────────────────────────────┤
│  services/       transport, retry, cache, logging   │  infrastructure
├─────────────────────────────────────────────────────┤
│  core/           sequence and protein computation   │  pure computation
│  models/         typed result dataclasses           │
└─────────────────────────────────────────────────────┘
```

Dependencies point downward only. `core` imports nothing from any layer above
it and nothing outside the standard library. This is enforced by a CI job
that installs the package with no extras and fails if `requests` or Biopython
can be imported.

---

## Why `core` has no dependencies

Three reasons, in order of importance.

**Testability.** Every calculation can be exercised offline, on any Python,
on any platform, in milliseconds. That is what makes a 432-test suite
practical and what lets the CI matrix cover 12 combinations without
installing a scientific stack.

**Installability.** The target user is a biologist who will type
`pip install` once. Requiring Biopython means requiring a build toolchain on
some machines, which is where installations fail.

**A stronger scientific claim.** If the package depended on Biopython, the
calculations would be Biopython's and their correctness assumed. Implementing
them independently and then *comparing* against Biopython
(`tests/test_biopython_parity.py`) demonstrates agreement instead of
assuming it. Where the two deliberately differ — ambiguity handling, Tm salt
correction defaults — the test documents why.

---

## Why results are dataclasses, not strings

Every v1.0 function returned a formatted string for display. That single
decision blocked batch processing, export, JSON output, testing of values
rather than wording, and reproducibility.

v2.0 functions return typed dataclasses (`models/results.py`) with
`.to_dict()`. Formatting happens once, at the presentation layer, in
`gui/formatting.py` — which is kept free of Tkinter imports precisely so the
displayed wording can be tested headless. The warnings that stop a
non-exact structural match being read as a confirmed one are part of that
tested wording.

---

## Why transports are pluggable

`services/transport.py` defines a transport interface with three
implementations: `requests`, `urllib`, and `FakeTransport`.

`FakeTransport` replays scripted responses, which is what makes failure-path
testing deterministic. Testing retry behaviour against a real service means
waiting for it to fail at the right moment; testing it against a script means
asserting that two 504s followed by a 200 produce exactly three attempts and
one recovery, on every machine, every time. The fault-injection benchmark is
the same machinery at a larger scale.

---

## Why the structure manager is separate from the providers

Each provider (`rcsb.py`, `alphafold.py`, `esmatlas.py`) knows one resource:
its endpoints, its quirks, what its 404 means. None of them knows about the
others, and none decides policy.

`manager.py` owns policy: which providers to try, in what order, when to
stop, what to do when one fails. It stops early **only** on an M4 exact
match. Anything less and it keeps going, because a related record from the
first resource is not a reason to skip a resource that might have the right
one.

Adding a provider means writing one class and adding its name to the fallback
order. See [`developer-guide.md`](developer-guide.md).

---

## Why validation is its own module

`structures/validation.py` contains no I/O and no knowledge of any specific
resource. It takes coordinate text, a query sequence, an accession and
settings, and returns a `ValidationReport`.

Keeping it separate means the M0–M4 classification is tested directly against
known inputs, independently of whether any service is reachable, and that
the same rules apply to every provider. A provider cannot exempt itself from
validation, which is the point.

---

## Threading

Network calls in the GUI run on a worker thread, with results returned to the
Tk main loop through a queue (`gui/widgets.py:BackgroundRunner`). Tk widgets
are only ever touched from the main thread.

Batch processing uses a thread pool. The work is I/O-bound — waiting on RCSB,
EBI and ESM Atlas — so threads are the right tool and the GIL is not a
constraint. Local-only batches are CPU-bound and run in the main thread,
which is why the offline throughput figure is single-threaded and honest.

Failures are isolated per record: one malformed sequence in a thousand is
reported in its own row and does not stop the run.

---

## Error handling

`services/errors.py` defines a typed hierarchy under `BioSeqInsightError`:
`TransportError`, `TimeoutError`, `NotFoundError`, `ServiceError`,
`ValidationError`, `CacheError`. Callers catch the category they can handle.

The rule throughout: **never substitute a different answer for a failed
one.** If the requested protein cannot be found, the software says so. It
does not return something similar and let the interface imply it is the same
thing. That behaviour is what the v1.0 audit found, and preventing it is the
reason this version exists.
