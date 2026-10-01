"""Resumable, rate-limited batches for large NCBI WGS project collections."""
import hashlib
import json
import os
import time
from pathlib import Path

from . import client

COLUMNS = ["batch", "projects", "accession", "description", "evalue", "bitscore", "identities_pct", "alignment_length", "subject_length", "subject_start", "subject_end", "query_alignment", "subject_alignment"]


def _write_json(path, value):
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(str(tmp), str(path))


def _write_tsv(path, rows):
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8", newline="") as handle:
        handle.write("\t".join(COLUMNS) + "\n")
        for row in rows:
            handle.write("\t".join(str(v).replace("\t", " ").replace("\n", " ") for v in row) + "\n")
    os.replace(str(tmp), str(path))


def run_batched_search(sequence, projects, output_dir, program="blastn", batch_size=250, expect=1e-5, hitlist_size=500, alignment_limit=500, max_batches=50, poll_interval=60, submission_interval=10, email=None, sleep=time.sleep):
    """Search each WGS project chunk and persist checkpoints so runs can resume."""
    normalized = list(dict.fromkeys(str(p).strip() for p in projects if str(p).strip()))
    if not normalized:
        raise ValueError("No WGS projects were supplied")
    if batch_size < 1 or max_batches < 1:
        raise ValueError("batch_size and max_batches must be positive")
    chunks = [normalized[i:i + batch_size] for i in range(0, len(normalized), batch_size)]
    if len(chunks) > max_batches:
        raise ValueError("Search needs {} batches; limit is {}. Increase batch_size or run a subset.".format(len(chunks), max_batches))

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    manifest_path = out / "manifest.json"
    manifest = {"sequence_sha256": hashlib.sha256(sequence.encode("utf-8")).hexdigest(), "projects": normalized, "program": program, "batch_size": batch_size, "expect": expect, "hitlist_size": hitlist_size, "alignment_limit": alignment_limit}
    if manifest_path.exists():
        if json.loads(manifest_path.read_text(encoding="utf-8")) != manifest:
            raise ValueError("Output directory belongs to a different search; choose another directory")
    else:
        if list(out.iterdir()):
            raise ValueError("Output directory is not empty and has no matching manifest")
        _write_json(manifest_path, manifest)

    for number, chunk in enumerate(chunks, 1):
        tsv = out / "batch_{:04d}.tsv".format(number)
        if tsv.exists():
            print("Batch {}/{} already complete; skipping".format(number, len(chunks)))
            continue
        checkpoint = out / "batch_{:04d}.rid.json".format(number)
        if checkpoint.exists():
            job = json.loads(checkpoint.read_text(encoding="utf-8"))
            rid = job["rid"]
            elapsed = time.time() - job["submitted_at"]
            if elapsed < poll_interval:
                sleep(poll_interval - elapsed)
        else:
            if number > 1 or submission_interval:
                sleep(submission_interval)
            params = {"EXPECT": expect, "HITLIST_SIZE": hitlist_size, "TOOL": "blast_api_wgs_batch"}
            if email:
                params["EMAIL"] = email
            rid = client.run_blast(sequence, programm=program, database=" ".join(chunk), wait=False, **params)
            _write_json(checkpoint, {"rid": rid, "submitted_at": time.time(), "batch": number, "projects": chunk})
            # NCBI asks clients to wait before the first status check, too.
            sleep(poll_interval)
        print("Batch {}/{} ({} projects), RID {}".format(number, len(chunks), len(chunk), rid))
        alignments = client.wait_for_blast_results(rid, poll_interval=poll_interval, alignment_limit=alignment_limit)
        rows = []
        for hit in alignments:
            start, end = hit.subj_range
            rows.append([number, len(chunk), hit.subj_id, hit.subj_name, hit.e_value, hit.score_bits, hit.identities, hit.align_len, hit.subj_len, start, end, hit.query_align, hit.sbjct_align])
        _write_tsv(tsv, rows)
        aggregate = []
        for batch_no in range(1, len(chunks) + 1):
            complete = out / "batch_{:04d}.tsv".format(batch_no)
            if complete.exists():
                lines = complete.read_text(encoding="utf-8").splitlines()
                aggregate.extend(line.split("\t") for line in lines[1:] if line)
        _write_tsv(out / "results.tsv", aggregate)
        print("Batch complete: {} hits; aggregate saved to {}".format(len(rows), out / "results.tsv"))
    return out / "results.tsv"
