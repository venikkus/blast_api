import pytest
from main import (
    Alignment,
    run_blast,
    parse_blast_text_output,
    extract_prefix_organism_pairs,
    filter_valid_wgs_ids,
    wait_for_blast_results
)
from pathlib import Path
import requests
import json
import logging
import pickle
import io


class TestAlignmentParsing:
    def test_query_alignments(self):
        blast_output = Path("tests/data/blast_output.txt").read_text(encoding="utf-8")
        actual = parse_blast_text_output(blast_output)
        expected = pickle.load(open("tests/data/parse_output.pkl", "rb"))
        assert actual[0].subj_id == expected[0].subj_id
        assert actual[99].subj_name == expected[99].subj_name

    def test_alignment_repr(self):
        aln = Alignment(
            subj_id="TEST123",
            subj_name="Test Organism",
            subj_len=1000,
            score_bits="99.9",
            score="200",
            e_value="1e-50",
            identities=90,
            match_len=90,
            align_len=100,
            query_align_chunks=[(1, "MKT", 3)],
            sbjct_align_chunks=[(1, "MKT", 3)],
        )
        rep = repr(aln)
        assert "TEST123" in rep
        assert "query_align='MKT'" in rep
        assert "subj_range=(1, 3)" in rep


class TestPrefixExtraction:
    def test_extract_prefix_organism_pairs(self):
        xml_text = Path("tests/data/prefixes_response.xml").read_text(encoding="utf-8")
        expected_lines = Path("tests/data/expected_prefixes.txt").read_text(encoding="utf-8").strip().splitlines()
        expected = [tuple(line.split("\t")) for line in expected_lines]
        result = extract_prefix_organism_pairs(xml_text)
        assert set(result) == set(expected)


class TestLogging:
    def test_logging_info_message(self):
        log_stream = io.StringIO()
        handler = logging.StreamHandler(log_stream)
        logger = logging.getLogger("test_logger")
        logger.setLevel(logging.INFO)
        logger.addHandler(handler)

        logger.info("Test message")
        handler.flush()

        log_contents = log_stream.getvalue()
        assert "Test message" in log_contents
        logger.removeHandler(handler)


class TestBlastIntegration:
    def test_run_blast_no_rid(self, monkeypatch):
        class MockResponse:
            def __init__(self):
                self.text = "RTOE = 10"
                self.status_code = 200

        def mock_post(url, data):
            return MockResponse()

        monkeypatch.setattr(requests, "post", mock_post)
        with pytest.raises(Exception, match="RID not found in response"):
            run_blast("FAKESEQ", wait=False)

    def test_wait_for_blast_results_failed(self, monkeypatch):
        def mock_get(url, params=None):
            class MockResponse:
                text = "Status=FAILED"
            return MockResponse()

        monkeypatch.setattr(requests, "get", mock_get)
        with pytest.raises(Exception, match="BLAST search failed"):
            wait_for_blast_results("FAKE_RID")


class TestWGSFiltering:
    def test_filter_valid_wgs_ids(self, monkeypatch):
        html_mock = Path("tests/data/db_species_table.html").read_text(encoding="utf-8")

        class MockResponse:
            def __init__(self, text, status_code=200):
                self.text = text
                self.status_code = status_code

        def mock_get(url, headers=None, params=None):
            return MockResponse(html_mock)

        monkeypatch.setattr(requests, "get", mock_get)

        prefix_list = ["ABJQ01", "AFNH01", "AFNH02", "JAHYSE01", "JAHYSF01", "JAQIFP01", "JBLWGC01"]
        expected = json.loads(Path("tests/data/expected_filtered_prefixes.json").read_text(encoding="utf-8"))
        result = filter_valid_wgs_ids(prefix_list)

        assert len(result) == len(expected)
        assert result == expected

    def test_filter_valid_wgs_ids_http_error(self, monkeypatch):
        def mock_get(url, headers=None, params=None):
            class MockResponse:
                status_code = 500
                text = "Internal Server Error"
            return MockResponse()

        monkeypatch.setattr(requests, "get", mock_get)
        with pytest.raises(Exception, match="Request failed with status code 500"):
            filter_valid_wgs_ids(["FAKE01", "FAKE02"])