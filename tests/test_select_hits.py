import csv
import hashlib
import json

import pytest

from blast_api.batching import COLUMNS
from blast_api.select_hits import read_query, select_hits, main


def hit(accession, evalue="1e-20", identity="95", alignment="AAAAAAAA", bitscore="100"):
    row = dict.fromkeys(COLUMNS, "")
    row.update(accession=accession, evalue=evalue, identities_pct=identity,
               query_alignment=alignment, bitscore=bitscore)
    return row


def batch(directory, number, rows):
    path = directory / "batch_{:04d}.tsv".format(number)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return path


def read(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_filters_and_ranks_each_batch_without_global_evalue_ranking(tmp_path):
    batch(tmp_path, 1, [hit("best", "1e-40"), hit("low_identity", "1e-80", "50"),
                        hit("short", "1e-80", alignment="AA"), hit("weak", "0.1")])
    batch(tmp_path, 2, [hit("second_batch", "1e-10")])
    out = tmp_path / "selected.tsv"
    stats = select_hits(tmp_path, out, 10, min_identity=90, min_coverage=80, top_per_batch=1)
    assert [r["accession"] for r in read(out)] == ["best", "second_batch"]
    assert stats == {"batches": 2, "rows": 5, "selected": 2}
    assert [r["rank_in_batch"] for r in read(out)] == ["1", "1"]


def test_gaps_and_distinct_accessions_and_extremely_small_evalues(tmp_path):
    batch(tmp_path, 1, [hit("x", "1e-400", alignment="AAA--AAAAA"),
                        hit("x", "1e-500"), hit("y", "e-450")])
    out = tmp_path / "selected.tsv"
    select_hits(tmp_path, out, 10, min_coverage=80, top_per_batch=0)
    rows = read(out)
    assert [r["evalue"] for r in rows] == ["1e-500", "e-450"]
    assert rows[0]["query_coverage_pct"] == "80.0000"
    assert rows[0]["coverage_method"] == "aligned_residues_estimate"


def test_coverage_sort_and_ties(tmp_path):
    batch(tmp_path, 1, [hit("short", "1e-100", alignment="AAA"), hit("long", "1e-10")])
    out = tmp_path / "selected.tsv"
    select_hits(tmp_path, out, 10, sort_by="coverage", top_per_batch=1)
    assert read(out)[0]["accession"] == "long"


def test_empty_batch_and_no_passing_hits_still_write_header(tmp_path):
    batch(tmp_path, 1, [])
    batch(tmp_path, 2, [hit("weak", "1")])
    out = tmp_path / "selected.tsv"
    assert select_hits(tmp_path, out, 10)["selected"] == 0
    assert read(out) == []
    assert "query_coverage_pct" in out.read_text()


@pytest.mark.parametrize("bad", ["NaN", "Infinity", "-1", "None"])
def test_bad_values_preserve_existing_output(tmp_path, bad):
    batch(tmp_path, 1, [hit("broken", bad)])
    out = tmp_path / "selected.tsv"
    out.write_text("previous results")
    with pytest.raises(ValueError, match="batch_0001.tsv:2"):
        select_hits(tmp_path, out, 10)
    assert out.read_text() == "previous results"


def test_coverage_over_100_is_preserved_for_review(tmp_path, capsys):
    batch(tmp_path, 1, [hit("overlap"), hit("valid", alignment="AAAA")])
    out = tmp_path / "selected.tsv"
    select_hits(tmp_path, out, 4)
    assert [row["accession"] for row in read(out)] == ["valid"]
    review = read(tmp_path / "selected.review.tsv")
    assert review[0]["accession"] == "overlap"
    assert review[0]["review_reason"] == "coverage_exceeds_100"
    assert "incomplete" in capsys.readouterr().err


def test_merged_hsps_below_100_percent_are_also_flagged(tmp_path):
    merged = hit("merged")
    merged["alignment_length"] = "4"
    batch(tmp_path, 1, [merged])
    out = tmp_path / "selected.tsv"
    select_hits(tmp_path, out, 10)
    assert read(out) == []
    assert read(tmp_path / "selected.review.tsv")[0]["review_reason"].startswith("alignment_length_mismatch")


def test_refuses_source_overwrite_and_missing_inputs(tmp_path):
    with pytest.raises(ValueError, match="No batch"):
        select_hits(tmp_path, tmp_path / "selected.tsv", 10)
    path = batch(tmp_path, 1, [hit("x")])
    with pytest.raises(ValueError, match="separate output"):
        select_hits(tmp_path, path, 10)


def test_query_and_manifest_validation(tmp_path):
    query = tmp_path / "query.fasta"
    query.write_text(">query\nACGT\n")
    assert read_query(query, tmp_path) == 4
    (tmp_path / "manifest.json").write_text(json.dumps({
        "sequence_sha256": hashlib.sha256(query.read_text().strip().encode()).hexdigest()}))
    assert read_query(query, tmp_path) == 4
    query.write_text(">query\nAAAA\n")
    with pytest.raises(ValueError, match="does not match"):
        read_query(query, tmp_path)
    query.write_text(">a\nAA\n>b\nTT\n")
    with pytest.raises(ValueError, match="Exactly one"):
        read_query(query, tmp_path)


def test_cli_accepts_percentages_and_rejects_invalid_thresholds(tmp_path):
    batch(tmp_path, 1, [hit("x")])
    out = tmp_path / "selected.tsv"
    assert main(["--input-dir", str(tmp_path), "--query-length", "10", "--output", str(out),
                 "--min-identity", "90", "--min-coverage", "80"]) == 0
    with pytest.raises(ValueError, match="Percentage thresholds"):
        select_hits(tmp_path, out, 10, min_coverage=101)


def test_defaults_keep_all_passing_accessions_in_api_and_cli(tmp_path):
    rows = [hit("hit{:02d}".format(i)) for i in range(15)]
    rows += [hit("low_identity", identity="40"), hit("short", alignment="AA")]
    batch(tmp_path, 1, rows)
    batch(tmp_path, 2, [hit("other_batch")])
    out = tmp_path / "selected.tsv"
    stats = select_hits(tmp_path, out, 10, min_identity=90, min_coverage=80)
    assert stats["selected"] == 16
    assert main(["--input-dir", str(tmp_path), "--query-length", "10",
                 "--output", str(out), "--min-identity", "90",
                 "--min-coverage", "80"]) == 0
    assert {row["accession"] for row in read(out)} == {
        "hit{:02d}".format(i) for i in range(15)
    } | {"other_batch"}
