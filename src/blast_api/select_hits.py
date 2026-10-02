"""Select hits independently per batch; legacy TSV coverage is an estimate."""
import argparse
import csv
import hashlib
import json
import os
import re
import sys
import tempfile
from decimal import Decimal, InvalidOperation
from pathlib import Path

EXTRA = ["source_batch", "rank_in_batch", "query_length", "query_coverage_pct", "coverage_method"]
REQUIRED = {"accession", "evalue", "bitscore", "identities_pct", "query_alignment"}


def number(value):
    text = str(value).strip()
    if text.lower().startswith("e"):
        text = "1" + text
    try:
        result = Decimal(text)
    except InvalidOperation:
        raise ValueError("Invalid number: {!r}".format(value))
    if not result.is_finite() or result < 0:
        raise ValueError("Expected a finite nonnegative number: {!r}".format(value))
    return result


def read_query(path, input_dir):
    text = Path(path).read_text(encoding="utf-8").strip()
    headers = sum(line.startswith(">") for line in text.splitlines())
    if headers > 1 or (headers and not text.startswith(">")):
        raise ValueError("Exactly one query sequence is supported")
    sequence = "".join(text.splitlines()[1:]) if headers else text
    sequence = "".join(sequence.split())
    if not re.fullmatch(r"[A-Za-z*]+", sequence):
        raise ValueError("Query must be one nonempty ungapped DNA/protein sequence")
    manifest = Path(input_dir) / "manifest.json"
    if manifest.exists():
        saved_hash = json.loads(manifest.read_text(encoding="utf-8")).get("sequence_sha256")
        submitted = text if headers else sequence
        if saved_hash and hashlib.sha256(submitted.encode("utf-8")).hexdigest() != saved_hash:
            raise ValueError("Query does not match manifest.json; use the original input FASTA")
    return len(sequence)


def select_hits(input_dir, output, query_length, max_evalue="1e-5", min_identity=0,
                min_coverage=0, top_per_batch=0, sort_by="evalue"):
    """Filter rows, rank within each file, and retain the best row per accession.

    Coverage counts nongap query residues, not the union of query coordinates:
    legacy batch tables lack coordinates and may contain merged HSPs.
    """
    max_evalue = number(max_evalue)
    min_identity, min_coverage = number(min_identity), number(min_coverage)
    if min_identity > 100 or min_coverage > 100:
        raise ValueError("Percentage thresholds must be between 0 and 100")
    if query_length < 1 or top_per_batch < 0:
        raise ValueError("query_length must be positive; top_per_batch must be nonnegative")
    if sort_by not in ("evalue", "coverage", "identity", "bitscore"):
        raise ValueError("Unknown sorting criterion")
    paths = sorted(p for p in Path(input_dir).glob("batch_*.tsv")
                   if re.fullmatch(r"batch_\d+\.tsv", p.name))
    if not paths:
        raise ValueError("No batch_NNNN.tsv files found in {}".format(input_dir))
    output = Path(output)
    if output.resolve() in {p.resolve() for p in paths} or output.name == "results.tsv":
        raise ValueError("Choose a separate output file, not a source batch or results.tsv")
    if re.fullmatch(r"batch_\d+\.tsv", output.name):
        raise ValueError("Output filename must not match batch_NNNN.tsv")
    selected, columns = [], None
    review = []
    review_path = output.with_name(output.stem + ".review.tsv")
    stats = {"batches": len(paths), "rows": 0, "selected": 0}
    for path in paths:
        candidates = []
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = reader.fieldnames or []
            if not REQUIRED.issubset(fields) or len(fields) != len(set(fields)):
                raise ValueError("{}: invalid header or missing required columns".format(path))
            if set(fields) & set(EXTRA):
                raise ValueError("{}: expected an original batch table".format(path))
            if columns is None:
                columns = fields
            elif fields != columns:
                raise ValueError("{}: batch headers differ".format(path))
            for row in reader:
                stats["rows"] += 1
                try:
                    if None in row or None in row.values() or not row["accession"].strip():
                        raise ValueError("Malformed row")
                    evalue, bitscore = number(row["evalue"]), number(row["bitscore"])
                    identity = number(row["identities_pct"])
                    if identity > 100:
                        raise ValueError("Identity exceeds 100%")
                    alignment = row["query_alignment"].strip()
                    if not re.fullmatch(r"[A-Za-z*\-]+", alignment):
                        raise ValueError("Invalid query_alignment")
                    covered = len(alignment.replace("-", ""))
                    coverage = Decimal(covered) * 100 / query_length
                    reason = None
                    if coverage > 100:
                        reason = "coverage_exceeds_100"
                    reported_length = row.get("alignment_length", "").strip()
                    if reported_length and number(reported_length) != len(alignment):
                        reason = "alignment_length_mismatch_possible_merged_HSPs"
                    if reason:
                        review.append(dict(row, source_batch=path.name,
                                           source_line=reader.line_num, review_reason=reason,
                                           query_length=query_length,
                                           aligned_query_residues=covered))
                        continue
                    if evalue > max_evalue or identity < min_identity or coverage < min_coverage:
                        continue
                    metrics = {"evalue": evalue, "coverage": -coverage,
                               "identity": -identity, "bitscore": -bitscore}
                    key = (metrics[sort_by], evalue, -coverage, -identity, -bitscore, row["accession"])
                    row.update(source_batch=path.name, query_length=query_length,
                               query_coverage_pct=format(coverage, ".4f"),
                               coverage_method="aligned_residues_estimate")
                    candidates.append((key, row))
                except ValueError as exc:
                    raise ValueError("{}:{}: {}".format(path, reader.line_num, exc))
        seen = set()
        rank = 0
        for _, row in sorted(candidates, key=lambda item: item[0]):
            if row["accession"] in seen:
                continue
            seen.add(row["accession"])
            rank += 1
            row["rank_in_batch"] = rank
            selected.append(row)
            if top_per_batch and rank >= top_per_batch:
                break
    # Validate every input before writing either output.
    write_table(review_path, columns + ["source_batch", "source_line", "review_reason",
                                       "query_length", "aligned_query_residues"], review)
    write_table(output, columns + EXTRA, selected)
    stats["selected"] = len(selected)
    if review:
        print("WARNING: {} ambiguous rows excluded from selection; preserved in {}. "
              "Selection is incomplete until these rows are reviewed.".format(
                  len(review), review_path), file=sys.stderr)
    return stats


def write_table(output, columns, rows):
    """Atomically replace one table, including a header for empty output."""
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", newline="",
                                         dir=str(output.parent), delete=False) as handle:
            temporary = Path(handle.name)
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(str(temporary), str(output))
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def main(argv=None):
    parser = argparse.ArgumentParser(description="Select the best hits from each WGS batch independently")
    parser.add_argument("--input-dir", required=True)
    query = parser.add_mutually_exclusive_group(required=True)
    query.add_argument("--query", help="Original single-sequence FASTA; checked against manifest when available")
    query.add_argument("--query-length", type=int, help="Query length in nucleotides (blastn) or amino acids (tblastn)")
    parser.add_argument("--output", required=True)
    parser.add_argument("--max-evalue", default="1e-5")
    parser.add_argument("--min-identity", type=float, default=0)
    parser.add_argument("--min-coverage", type=float, default=0)
    parser.add_argument("--top-per-batch", type=int, default=0, help="Optional limit per batch; default 0 = all passing accessions")
    parser.add_argument("--sort-by", choices=("evalue", "coverage", "identity", "bitscore"), default="evalue")
    args = parser.parse_args(argv)
    try:
        length = read_query(args.query, args.input_dir) if args.query else args.query_length
        print("Coverage is an estimate from nongap query residues; merged overlapping HSPs can overestimate it.", file=sys.stderr)
        stats = select_hits(args.input_dir, args.output, length, args.max_evalue,
                            args.min_identity, args.min_coverage, args.top_per_batch, args.sort_by)
    except (OSError, ValueError, csv.Error) as exc:
        parser.exit(2, "Error: {}\n".format(exc))
    print("Batches: {batches}; rows read: {rows}; hits selected: {selected}".format(**stats))
    print("Saved: {}".format(args.output))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
