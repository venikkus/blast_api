"""Command line entry point for batched WGS searches."""
import argparse
import re
import sys
import time
from pathlib import Path
from .batching import run_batched_search
from .client import get_wgs_projects, resolve_taxon_id


def _read_fasta(path):
    text = Path(path).read_text(encoding="utf-8").strip()
    if not text:
        raise ValueError("Query FASTA is empty")
    if text.startswith(">"):
        sequence = "".join(line.strip() for line in text.splitlines() if not line.startswith(">"))
        if not sequence:
            raise ValueError("FASTA contains no sequence")
        return text
    sequence = re.sub(r"\s+", "", text)
    if not sequence:
        raise ValueError("FASTA contains no sequence")
    return sequence


def make_parser():
    parser = argparse.ArgumentParser(prog="blast-wgs", description="Search selected NCBI WGS projects in resumable batches.")
    parser.add_argument("--query", required=True, help="Input FASTA file")
    parser.add_argument("--organism", action="append", default=[], help="Taxon name; repeatable")
    parser.add_argument("--taxid", action="append", default=[], help="NCBI TaxID; repeatable")
    parser.add_argument("--exclude-taxid", action="append", default=[], help="TaxID to exclude")
    parser.add_argument("--program", choices=("blastn", "tblastn"), default="blastn")
    parser.add_argument("--batch-size", type=int, default=250)
    parser.add_argument("--max-batches", type=int, default=50)
    parser.add_argument("--expect", type=float, default=1e-5)
    parser.add_argument("--hitlist-size", type=int, default=500)
    parser.add_argument("--alignment-limit", type=int, default=500)
    parser.add_argument("--output-dir", default="results/wgs_run")
    parser.add_argument("--email", help="Optional contact email for NCBI")
    return parser


def main(argv=None):
    args = make_parser().parse_args(argv)
    try:
        if not args.organism and not args.taxid:
            raise ValueError("Pass at least one --organism or --taxid")
        taxids = list(args.taxid)
        for index, organism in enumerate(args.organism):
            taxids.append(organism.strip() if organism.strip().isdigit() else resolve_taxon_id(organism))
            if index + 1 < len(args.organism):
                time.sleep(10)
        time.sleep(10)
        projects = get_wgs_projects(taxids, exclude_taxids=args.exclude_taxid or None)
        print("NCBI returned {} WGS projects".format(len(projects)))
        if not projects:
            raise ValueError("No WGS projects found for selected taxa")
        print("They will be searched in batches of {}; {} batch(es) maximum per run.".format(args.batch_size, args.max_batches))
        print("Batch E-values are calculated independently and should not be compared directly.")
        result = run_batched_search(_read_fasta(args.query), projects, args.output_dir, program=args.program, batch_size=args.batch_size, expect=args.expect, hitlist_size=args.hitlist_size, alignment_limit=args.alignment_limit, max_batches=args.max_batches, email=args.email)
        print("Finished: {}".format(result))
        return 0
    except (OSError, ValueError, LookupError) as error:
        print("Error: {}".format(error), file=sys.stderr)
        return 2
