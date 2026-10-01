import json
from unittest.mock import Mock, patch

import pytest

from blast_api.batching import run_batched_search


class Hit:
    subj_id = "ABC123.1"
    subj_name = "example"
    e_value = "1e-20"
    score_bits = "90"
    identities = 90
    align_len = 100
    subj_len = 120
    subj_range = (1, 100)
    query_align = "A" * 100
    sbjct_align = "A" * 100


def test_run_batched_search_saves_rids_and_aggregates(tmp_path):
    with patch("blast_api.batching.client.run_blast", side_effect=["RID1", "RID2"]) as submit, \
            patch("blast_api.batching.client.wait_for_blast_results", return_value=[Hit()]) as wait:
        output = run_batched_search("ACGT", ["WGS_VDB://AAAA01", "WGS_VDB://BBBB01", "WGS_VDB://CCCC01"], tmp_path, batch_size=2, submission_interval=0, sleep=lambda _seconds: None)

    assert output.exists()
    assert len(output.read_text().splitlines()) == 3
    assert submit.call_count == 2
    assert wait.call_count == 2
    assert json.loads((tmp_path / "batch_0001.rid.json").read_text())["rid"] == "RID1"


def test_run_batched_search_resumes_checkpoint_without_resubmitting(tmp_path):
    with patch("blast_api.batching.client.run_blast", return_value="RID1"), \
            patch("blast_api.batching.client.wait_for_blast_results", return_value=[]):
        run_batched_search("ACGT", ["WGS_VDB://AAAA01"], tmp_path, submission_interval=0, poll_interval=0)
    with patch("blast_api.batching.client.run_blast") as submit, \
            patch("blast_api.batching.client.wait_for_blast_results", return_value=[]):
        run_batched_search("ACGT", ["WGS_VDB://AAAA01"], tmp_path, submission_interval=0, poll_interval=0)
    submit.assert_not_called()


def test_run_batched_search_rejects_changed_search_parameters(tmp_path):
    with patch("blast_api.batching.client.run_blast", return_value="RID1"), \
            patch("blast_api.batching.client.wait_for_blast_results", return_value=[]):
        run_batched_search("ACGT", ["WGS_VDB://AAAA01"], tmp_path, submission_interval=0, poll_interval=0)
    with pytest.raises(ValueError, match="different search"):
        run_batched_search("TGCA", ["WGS_VDB://AAAA01"], tmp_path, submission_interval=0)


def test_run_batched_search_checks_max_batches_before_submitting(tmp_path):
    with patch("blast_api.batching.client.run_blast") as submit:
        with pytest.raises(ValueError, match="needs 3 batches"):
            run_batched_search("ACGT", ["A", "B", "C"], tmp_path, batch_size=1, max_batches=2)
    submit.assert_not_called()
