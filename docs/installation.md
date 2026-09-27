# Installation

BioSeqInsight needs **Python 3.10 or newer** and nothing else. No compiler,
no scientific stack, no system libraries.

---

## Quick install

```bash
pip install bioseqinsight
bioseqinsight version
```

## From source

```bash
git clone https://github.com/niloy-999/BioSeqInsight.git
cd BioSeqInsight
python -m venv .venv
source .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -e ".[dev]"
python -m pytest -q
```

---

## Optional extras

| Extra | Installs | Why you might want it |
|---|---|---|
| `http` | `requests` | Connection pooling and better proxy handling. The stdlib transport is used automatically when absent, and the software behaves identically either way. |
| `validation` | `biopython` | Enables `tests/test_biopython_parity.py`, which cross-validates every calculation against Biopython. Needed to reproduce the parity results in the paper, not to run the software. |
| `dev` | pytest, ruff, black, mypy, build | Development and CI. |

```bash
pip install "bioseqinsight[http]"
pip install "bioseqinsight[http,validation]"
```

---

## The graphical interface

The GUI needs Tkinter, which ships with most Python builds.

| Platform | If Tkinter is missing |
|---|---|
| Debian / Ubuntu | `sudo apt install python3-tk` |
| Fedora / RHEL | `sudo dnf install python3-tkinter` |
| Arch | `sudo pacman -S tk` |
| macOS (python.org installer) | included |
| macOS (Homebrew) | `brew install python-tk` |
| Windows | included; re-run the installer and tick "tcl/tk and IDLE" if not |

Check:

```bash
python -c "import tkinter; print(tkinter.TkVersion)"
bioseqinsight version        # reports whether Tkinter is available
```

**Tkinter is not required.** Every capability is available from the command
line, so headless servers, containers and CI are fully supported. Launching
the GUI without Tkinter prints an explanation and the CLI alternative rather
than a traceback.

---

## Verifying the installation

```bash
bioseqinsight version
bioseqinsight --offline dna ATGGCTAGCTAGGCTTACCGATTGGCATTAA --min-orf 3
bioseqinsight --offline protein MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEG
python -m pytest -q          # from a source checkout
```

The test suite runs entirely offline. If it passes, the installation is
sound.

To check network access to the structural resources:

```bash
bioseqinsight structure P0CG48
```

Expect an M4 exact match to PDB entry 1UBQ.

---

## Platform notes

**Windows.** Use `py -m pip install bioseqinsight` if `pip` is not on PATH.
Console scripts land in `%LOCALAPPDATA%\Programs\Python\Python3xx\Scripts`;
add it to PATH or use `py -m bioseqinsight`.

**macOS.** On Apple Silicon everything is pure Python, so there is no
architecture-specific build step. If `bioseqinsight` is not found after
install, `python3 -m bioseqinsight` always works.

**Linux without a display.** Use the CLI. `bioseqinsight-gui` will exit with
an explanatory message rather than crashing.

**Behind a proxy.** The stdlib transport honours `HTTPS_PROXY` and
`NO_PROXY`. Installing the `http` extra gives `requests`, which handles
corporate proxy configurations more flexibly.

---

## Configuration file location

| Platform | Default |
|---|---|
| Linux / macOS | `~/.bioseqinsight/settings.json` |
| Windows | `%USERPROFILE%\.bioseqinsight\settings.json` |

Create one with the current effective settings:

```bash
bioseqinsight config --save ~/.bioseqinsight/settings.json
```

Downloads, cache and logs default to sibling directories under
`~/.bioseqinsight/` and are all configurable.

---

## Troubleshooting

**`command not found: bioseqinsight`** — the scripts directory is not on
PATH. Use `python -m bioseqinsight` instead, which always works.

**`ModuleNotFoundError: No module named 'bioseqinsight'`** — you are in a
different environment from the one you installed into. Check with
`python -c "import sys; print(sys.executable)"`.

**`No module named 'tkinter'`** — see the table above. The CLI is unaffected.

**Network timeouts** — public services do go down. Raise the retry budget
with `BIOSEQINSIGHT_MAX_ATTEMPTS=5`, or work from the cache. `--offline`
disables all outbound requests.

**A proxy intercepts HTTPS with a private CA** — set `SSL_CERT_FILE` to your
CA bundle. The stdlib transport refuses plain HTTP by design, so a proxy that
downgrades the connection will fail rather than transmit in the clear.

---

## Uninstalling

```bash
pip uninstall bioseqinsight
rm -rf ~/.bioseqinsight      # cache, logs, settings and downloads
```

Project directories you created are ordinary folders and are left alone.
