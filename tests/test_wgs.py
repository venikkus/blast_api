import json
from unittest.mock import Mock, patch

import pytest

from blast_api import client


def test_get_wgs_projects_parses_and_deduplicates_ncbi_response():
    response = Mock()
    response.text = "WGS_VDB://JPYT01\nWGS_VDB://JPYT01 WGS_VDB://JZAA01\n"

    with patch("blast_api.client.requests.post", return_value=response) as post:
        projects = client.get_wgs_projects(["9606", 10090], exclude_taxids=[10091])

    assert projects == ["WGS_VDB://JPYT01", "WGS_VDB://JZAA01"]
    post.assert_called_once_with(
        "https://www.ncbi.nlm.nih.gov/blast/BDB2EZ/taxid2wgs.cgi",
        params={"INCLUDE_TAXIDS": "9606,10090", "EXCLUDE_TAXIDS": "10091"},
        timeout=60,
    )


def test_get_wgs_projects_rejects_non_numeric_taxids_without_request():
    with patch("blast_api.client.requests.post") as post:
        with pytest.raises(ValueError, match="numeric NCBI TaxIDs"):
            client.get_wgs_projects(["9606", "Homo sapiens"])

    post.assert_not_called()


def test_resolve_taxon_id_returns_unique_ncbi_match():
    response = Mock()
    response.content = b"<eSearchResult><IdList><Id>9606</Id></IdList></eSearchResult>"

    with patch("blast_api.client.requests.get", return_value=response) as get:
        taxid = client.resolve_taxon_id("Homo sapiens")

    assert taxid == "9606"
    assert get.call_args.kwargs["params"] == {
        "db": "taxonomy",
        "term": '"Homo sapiens"[name]',
        "retmode": "xml",
        "retmax": 20,
    }


def test_resolve_taxon_id_reports_ambiguous_names():
    search_response = Mock()
    search_response.content = (
        b"<eSearchResult><IdList><Id>111</Id><Id>222</Id></IdList></eSearchResult>"
    )
    summary_response = Mock()
    summary_response.json.return_value = {
        "result": {
            "111": {"taxname": "Example alpha", "rank": "species"},
            "222": {"taxname": "Example beta", "rank": "genus"},
        }
    }

    with patch("blast_api.client.requests.get", side_effect=[search_response, summary_response]):
        with pytest.raises(ValueError, match="ambiguous") as error:
            client.resolve_taxon_id("Example")

    assert "Example alpha" in str(error.value)
    assert "TaxID 222" in str(error.value)


def test_run_blast_resolves_taxon_name_before_submitting_wgs_search():
    blast_response = Mock()
    blast_response.status_code = 200
    blast_response.text = "RID = TEST123\nRTOE = 1\n"

    with patch("blast_api.client.resolve_taxon_id", return_value="9606") as resolve, \
            patch("blast_api.client.get_wgs_projects", return_value=["WGS_VDB://JPYT01"]), \
            patch("blast_api.client.requests.post", return_value=blast_response) as post:
        rid = client.run_blast(
            "ACGT", programm="blastn", database="wgs", taxon="Homo sapiens",
            wait=False,
        )

    assert rid == "TEST123"
    resolve.assert_called_once_with("Homo sapiens")
    assert post.call_args.kwargs["data"] == {
        "CMD": "Put",
        "PROGRAM": "blastn",
        "DATABASE": "WGS_VDB://JPYT01",
        "QUERY": "ACGT",
        "ENTREZ_QUERY": None,
    }


def test_parser_keeps_hsps_separate_with_their_own_metrics():
    text = '''>ABC123 example subject
Length=1000
 Score = 100 bits (200), Expect = 1e-30
 Identities = 4/4 (100%)
Query  1  ACGT  4
          ||||
Sbjct  10  ACGT  13

 Score = 50 bits (90), Expect = 1e-10
 Identities = 3/4 (75%)
Query  1  ACGT  4
          |||
Sbjct  30  ACGA  33
'''
    hits = client.parse_blast_text_output(text)
    assert len(hits) == 2
    assert [h.query_align for h in hits] == ["ACGT", "ACGT"]
    assert [h.e_value for h in hits] == ["1e-30", "1e-10"]
    assert [h.identities for h in hits] == [100, 75]
    assert [h.subj_len for h in hits] == [1000, 1000]
    assert [h.subj_range for h in hits] == [(10, 13), (30, 33)]
