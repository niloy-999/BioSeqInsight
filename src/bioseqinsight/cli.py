"""Command-line interface.

Everything the GUI can do is reachable from the command line. That is not a
convenience feature: it is what makes the evaluation in the manuscript
reproducible. A reviewer can clone the repository, run one command and
regenerate a benchmark table without a display server, a mouse or a
screenshot.

Subcommands
-----------
``dna``        analyse a nucleotide sequence or FASTA file
``protein``    analyse an amino-acid sequence or FASTA file
``structure``  retrieve and validate structures for a sequence/accession/PDB id
``batch``      run a whole FASTA file and write CSV/TSV/JSON/HTML
``project``    create, inspect, export and import project directories
``cache``      inspect or clear the response cache
``config``     show or write the effective settings
``version``    print version and environment information
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import __version__
from .config.settings import Settings
from .core.alphabet import SequenceError, parse_fasta
from .core.protein import ProteinAnalyzer
from .core.sequence import SequenceAnalyzer
from .io.export import (
    to_json,
    write_csv,
    write_html_report,
    write_json,
    write_tsv,
)
from .services.logging_setup import configure_logging
from .workflows.batch import BatchOptions, BatchRunner
from .workflows.project import Project, ProjectError, list_projects

EXIT_OK = 0
EXIT_USER_ERROR = 2
EXIT_RUNTIME_ERROR = 3


# --------------------------------------------------------------------------
# argument parsing
# --------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="bioseqinsight",
        description=(
            "Integrated DNA sequence analysis and validated protein "
            "structure-resource retrieval."
        ),
    )
    parser.add_argument("--version", action="version", version=f"BioSeqInsight {__version__}")
    parser.add_argument("--config", help="path to a JSON settings file")
    parser.add_argument("--offline", action="store_true", help="never contact external services")
    parser.add_argument("--no-cache", action="store_true", help="bypass the response cache")
    parser.add_argument("--refresh", action="store_true", help="ignore cached entries and refetch")
    parser.add_argument("--download-dir", help="where retrieved coordinates are written")
    parser.add_argument("--log-level", default=None, choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    parser.add_argument("--quiet", action="store_true", help="suppress log output on stderr")

    sub = parser.add_subparsers(dest="command", required=True)

    dna = sub.add_parser("dna", help="analyse a nucleotide sequence")
    _add_input_args(dna)
    dna.add_argument("--motif", help="IUPAC motif to search for")
    dna.add_argument("--min-orf", type=int, default=30, help="minimum ORF length in amino acids")
    dna.add_argument("--table", type=int, default=1, help="NCBI genetic code table")
    dna.add_argument("--json", action="store_true", help="emit JSON instead of text")
    dna.add_argument("--out", help="write the result to this file")

    protein = sub.add_parser("protein", help="analyse an amino-acid sequence")
    _add_input_args(protein)
    protein.add_argument("--json", action="store_true")
    protein.add_argument("--sketch", action="store_true", help="include the propensity sketch")
    protein.add_argument("--out")

    structure = sub.add_parser("structure", help="retrieve and validate a structure")
    structure.add_argument("query", help="protein sequence, UniProt accession or PDB identifier")
    structure.add_argument(
        "--providers",
        help="comma-separated subset of rcsb,alphafold,esmatlas (default: all, in settings order)",
    )
    structure.add_argument("--all", action="store_true", help="query every provider, do not stop at M4")
    structure.add_argument("--json", action="store_true")
    structure.add_argument("--out")

    batch = sub.add_parser("batch", help="analyse a FASTA file and export a table")
    batch.add_argument("fasta", help="input FASTA file")
    batch.add_argument("--structures", action="store_true", help="also retrieve structures")
    batch.add_argument("--providers")
    batch.add_argument("--workers", type=int, default=4)
    batch.add_argument("--min-orf", type=int, default=30)
    batch.add_argument("--alphabet", choices=["dna", "protein"], help="skip autodetection")
    batch.add_argument("--out-prefix", default="bioseqinsight_batch")
    batch.add_argument(
        "--formats",
        default="csv,json",
        help="comma-separated subset of csv,tsv,json,html",
    )
    batch.add_argument("--project", help="write results into this project directory")

    project = sub.add_parser("project", help="manage project directories")
    project_sub = project.add_subparsers(dest="project_command", required=True)
    create = project_sub.add_parser("create")
    create.add_argument("path")
    create.add_argument("--name", default="")
    create.add_argument("--description", default="")
    create.add_argument("--fasta", help="import this FASTA file immediately")
    info = project_sub.add_parser("info")
    info.add_argument("path")
    export = project_sub.add_parser("export")
    export.add_argument("path")
    export.add_argument("--out")
    imp = project_sub.add_parser("import")
    imp.add_argument("archive")
    imp.add_argument("destination")
    listing = project_sub.add_parser("list")
    listing.add_argument("root", nargs="?", default=None)

    cache = sub.add_parser("cache", help="inspect or clear the response cache")
    cache_sub = cache.add_subparsers(dest="cache_command", required=True)
    cache_sub.add_parser("stats")
    cache_list = cache_sub.add_parser("list")
    cache_list.add_argument("--namespace")
    cache_clear = cache_sub.add_parser("clear")
    cache_clear.add_argument("--namespace")

    config = sub.add_parser("config", help="show or save the effective settings")
    config.add_argument("--save", help="write the effective settings to this path")

    sub.add_parser("version", help="print version and environment details")
    return parser


def _add_input_args(parser: argparse.ArgumentParser) -> None:
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("sequence", nargs="?", help="sequence text, or '-' to read stdin")
    group.add_argument("--file", help="FASTA file (the first record is used)")


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------


def _settings_from_args(args: argparse.Namespace) -> Settings:
    overrides: dict = {}
    if args.offline:
        overrides["offline"] = True
    if args.no_cache:
        overrides["cache_enabled"] = False
    if args.download_dir:
        overrides["download_dir"] = args.download_dir
    if args.log_level:
        overrides["log_level"] = args.log_level
    return Settings.load(path=args.config, **overrides)


def _read_input(args: argparse.Namespace) -> tuple[str, str]:
    """Return ``(identifier, sequence_text)`` from the CLI's input options."""
    if getattr(args, "file", None):
        path = Path(args.file)
        if not path.is_file():
            raise SequenceError(f"No such file: {path}")
        records = parse_fasta(path.read_text(encoding="utf-8"))
        if not records:
            raise SequenceError(f"{path} contains no sequence records.")
        if len(records) > 1:
            print(
                f"note: {path.name} holds {len(records)} records; analysing the first. "
                "Use 'bioseqinsight batch' for all of them.",
                file=sys.stderr,
            )
        return records[0].identifier, records[0].sequence
    text = args.sequence
    if text == "-":
        text = sys.stdin.read()
    records = parse_fasta(text)
    if records:
        return records[0].identifier, records[0].sequence
    return "sequence", text or ""


def _emit(text: str, out: str | None) -> None:
    if out:
        Path(out).parent.mkdir(parents=True, exist_ok=True)
        Path(out).write_text(text, encoding="utf-8")
        print(f"Written to {out}")
    else:
        print(text)


# --------------------------------------------------------------------------
# commands
# --------------------------------------------------------------------------


def cmd_dna(args: argparse.Namespace, settings: Settings) -> int:
    identifier, text = _read_input(args)
    analyzer = SequenceAnalyzer(text, identifier=identifier, table=args.table)
    result = analyzer.analyze(motif=args.motif, min_orf_aa=args.min_orf)
    if args.json:
        _emit(to_json(result), args.out)
        return EXIT_OK

    lines = [
        f"Sequence: {result.identifier}",
        f"Length: {result.length} nt",
        f"GC content: {result.gc_percent:.2f}%",
        "Composition: "
        + ", ".join(f"{base}={count}" for base, count in result.composition.items()),
    ]
    if result.tm:
        lines.append(f"Tm: {result.tm.tm_c:.2f} C ({result.tm.method})")
        if result.tm.note:
            lines.append(f"  note: {result.tm.note}")
    lines.append(f"ORFs (>= {args.min_orf} aa): {len(result.orfs)}")
    for orf in result.orfs[:5]:
        lines.append(
            f"  strand {orf.strand} frame {orf.frame} nt {orf.start}-{orf.end} "
            f"({orf.aa_length} aa){' [no stop]' if orf.unterminated else ''}"
        )
    if args.motif:
        lines.append(f"Motif {args.motif}: {len(result.motifs)} hit(s)")
        for hit in result.motifs[:10]:
            lines.append(f"  {hit.strand} strand at {hit.start}-{hit.end}")
    for warning in result.warnings:
        lines.append(f"warning: {warning}")
    _emit("\n".join(lines), args.out)
    return EXIT_OK


def cmd_protein(args: argparse.Namespace, settings: Settings) -> int:
    identifier, text = _read_input(args)
    result = ProteinAnalyzer(text, identifier=identifier).analyze(include_sketch=args.sketch)
    if args.json:
        _emit(to_json(result), args.out)
        return EXIT_OK
    lines = [
        f"Protein: {result.identifier}",
        f"Length: {result.length} aa",
        f"Molecular weight: {result.molecular_weight:.2f} Da ({result.molecular_weight / 1000:.3f} kDa)",
        f"GRAVY: {result.gravy:.4f}",
        f"Isoelectric point: {result.isoelectric_point:.2f}",
        f"Aromaticity: {result.aromaticity:.4f}",
        f"Aliphatic index: {result.aliphatic_index:.2f}",
        f"Extinction (280 nm): {result.extinction_coefficient_reduced} reduced / "
        f"{result.extinction_coefficient_cystines} with cystines",
    ]
    if result.secondary_structure:
        sketch = result.secondary_structure
        lines += [
            "",
            f"Propensity sketch ({sketch.method}):",
            f"  helix {sketch.helix_percent:.1f}%  sheet {sketch.sheet_percent:.1f}%  "
            f"coil {sketch.coil_percent:.1f}%",
            f"  {sketch.visual}",
        ]
    for warning in result.warnings:
        lines.append(f"warning: {warning}")
    _emit("\n".join(lines), args.out)
    return EXIT_OK


def cmd_structure(args: argparse.Namespace, settings: Settings) -> int:
    from .structures.manager import StructureManager  # imported late: pulls in networking

    manager = StructureManager(settings)
    providers = args.providers.split(",") if args.providers else None
    outcome = manager.retrieve(
        args.query,
        providers=providers,
        stop_on_exact=not args.all,
        refresh=args.refresh,
    )
    if args.json:
        _emit(to_json(outcome), args.out)
        return EXIT_OK if outcome.best else EXIT_RUNTIME_ERROR

    lines = [
        f"Query: {outcome.query.raw[:60]}{'...' if len(outcome.query.raw) > 60 else ''}",
        f"Interpreted as: {outcome.query.query_type}",
    ]
    if outcome.query.accession:
        lines.append(f"Accession: {outcome.query.accession}")
    lines.append("")
    for result in outcome.results:
        lines.append(f"[{result.source}] {result.status.value}")
        if result.identifier:
            lines.append(f"  identifier: {result.identifier}")
        if result.download_path:
            lines.append(f"  saved: {result.download_path}")
        if result.confidence is not None:
            lines.append(f"  {result.confidence_kind}: {result.confidence:.2f}")
        if result.attempts:
            lines.append(
                f"  attempts: {len(result.attempts)}"
                + (" (recovered after failure)" if result.recovered else "")
            )
        if result.from_cache:
            lines.append("  served from cache")
        validation = result.validation
        if validation:
            lines.append(f"  mapping: {validation.mapping_level.value} {validation.mapping_level.label}")
            if validation.identity_percent is not None:
                lines.append(
                    f"  identity {validation.identity_percent:.1f}%  "
                    f"coverage {validation.coverage_percent:.1f}%"
                )
            if validation.message:
                lines.append(f"  {validation.message}")
        if result.error:
            lines.append(f"  error: {result.error}")
        lines.append("")

    best = outcome.best
    lines.append(
        f"Best result: {best.source} {best.identifier} ({outcome.best_mapping.value} "
        f"{outcome.best_mapping.label})"
        if best
        else "Best result: none"
    )
    for warning in outcome.warnings:
        lines.append(f"warning: {warning}")
    _emit("\n".join(lines), args.out)
    return EXIT_OK if best else EXIT_RUNTIME_ERROR


def cmd_batch(args: argparse.Namespace, settings: Settings) -> int:
    path = Path(args.fasta)
    if not path.is_file():
        print(f"No such FASTA file: {path}", file=sys.stderr)
        return EXIT_USER_ERROR

    runner = BatchRunner(settings)
    records = runner.records_from_fasta_file(str(path))
    if not records:
        print(f"{path} contains no sequence records.", file=sys.stderr)
        return EXIT_USER_ERROR

    options = BatchOptions(
        include_structures=args.structures,
        providers=args.providers.split(",") if args.providers else None,
        min_orf_aa=args.min_orf,
        workers=args.workers,
        refresh_cache=args.refresh,
        force_alphabet=args.alphabet,
    )

    def progress(done: int, total: int, message: str) -> None:
        print(f"\r[{done}/{total}] {message[:60]:<60}", end="", file=sys.stderr, flush=True)

    rows, summary = runner.run(records, options, progress=None if args.quiet else progress)
    if not args.quiet:
        print(file=sys.stderr)

    formats = {f.strip() for f in args.formats.split(",") if f.strip()}
    written: dict[str, str] = {}
    if args.project:
        project = Project.open(args.project) if Path(args.project, "project.json").is_file() \
            else Project.create(args.project, settings=settings)
        project.add_sequences(records, filename=path.name)
        written = project.save_results(rows, basename=args.out_prefix, formats=sorted(formats))
    else:
        if "csv" in formats:
            written["csv"] = write_csv(rows, f"{args.out_prefix}.csv")
        if "tsv" in formats:
            written["tsv"] = write_tsv(rows, f"{args.out_prefix}.tsv")
        if "json" in formats:
            written["json"] = write_json(
                {
                    "software_version": __version__,
                    "settings": settings.to_dict(),
                    "summary": summary.to_dict(),
                    "rows": [row.to_dict() for row in rows],
                },
                f"{args.out_prefix}.json",
            )
        if "html" in formats:
            written["html"] = write_html_report(rows, f"{args.out_prefix}.html")

    print(json.dumps(summary.to_dict(), indent=2))
    for kind, target in written.items():
        print(f"{kind}: {target}")
    return EXIT_OK


def cmd_project(args: argparse.Namespace, settings: Settings) -> int:
    try:
        if args.project_command == "create":
            project = Project.create(
                args.path, name=args.name, description=args.description, settings=settings
            )
            if args.fasta:
                project.import_fasta(args.fasta)
            print(json.dumps(project.summary(), indent=2))
        elif args.project_command == "info":
            print(json.dumps(Project.open(args.path).summary(), indent=2))
        elif args.project_command == "export":
            print(Project.open(args.path).export_zip(args.out))
        elif args.project_command == "import":
            project = Project.import_zip(args.archive, args.destination)
            print(json.dumps(project.summary(), indent=2))
        elif args.project_command == "list":
            root = args.root or settings.projects_dir
            print(json.dumps(list_projects(root), indent=2))
    except ProjectError as exc:
        print(str(exc), file=sys.stderr)
        return EXIT_USER_ERROR
    return EXIT_OK


def cmd_cache(args: argparse.Namespace, settings: Settings) -> int:
    from .services.cache import ResponseCache

    cache = ResponseCache(settings.cache_dir, ttl_s=settings.cache_ttl_s)
    if args.cache_command == "stats":
        print(json.dumps(cache.stats(), indent=2))
    elif args.cache_command == "list":
        entries = cache.entries(args.namespace)
        print(
            json.dumps(
                [
                    {
                        "namespace": e.namespace,
                        "query": e.query[:60],
                        "source": e.source,
                        "created_at": e.created_at,
                        "size_bytes": e.size_bytes,
                        "url": e.url,
                    }
                    for e in entries
                ],
                indent=2,
            )
        )
    elif args.cache_command == "clear":
        print(f"Removed {cache.clear(args.namespace)} cached payload(s).")
    return EXIT_OK


def cmd_config(args: argparse.Namespace, settings: Settings) -> int:
    if args.save:
        print(f"Settings written to {settings.save(args.save)}")
        return EXIT_OK
    print(json.dumps({"source": settings.source, **settings.to_dict()}, indent=2, sort_keys=True))
    return EXIT_OK


def cmd_version(args: argparse.Namespace, settings: Settings) -> int:
    import platform

    optional = {}
    for module in ("Bio", "requests"):
        try:
            imported = __import__(module)
            optional[module] = getattr(imported, "__version__", "installed")
        except ImportError:
            optional[module] = "not installed"
    try:
        import tkinter  # noqa: F401

        optional["tkinter"] = "available"
    except ImportError:
        optional["tkinter"] = "not available (GUI disabled)"

    print(
        json.dumps(
            {
                "bioseqinsight": __version__,
                "python": platform.python_version(),
                "platform": platform.platform(),
                "optional_dependencies": optional,
                "settings_source": settings.source,
            },
            indent=2,
        )
    )
    return EXIT_OK


COMMANDS = {
    "dna": cmd_dna,
    "protein": cmd_protein,
    "structure": cmd_structure,
    "batch": cmd_batch,
    "project": cmd_project,
    "cache": cmd_cache,
    "config": cmd_config,
    "version": cmd_version,
}


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        settings = _settings_from_args(args)
    except ValueError as exc:
        print(f"Configuration error: {exc}", file=sys.stderr)
        return EXIT_USER_ERROR

    configure_logging(
        log_file=settings.log_file,
        level=settings.log_level,
        json_format=settings.log_json,
        console=not args.quiet,
    )
    try:
        return COMMANDS[args.command](args, settings)
    except SequenceError as exc:
        print(f"Input error: {exc}", file=sys.stderr)
        return EXIT_USER_ERROR
    except KeyboardInterrupt:  # pragma: no cover - interactive
        print("\nInterrupted.", file=sys.stderr)
        return EXIT_RUNTIME_ERROR
    except Exception as exc:  # pragma: no cover - last resort
        print(f"{type(exc).__name__}: {exc}", file=sys.stderr)
        return EXIT_RUNTIME_ERROR


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
