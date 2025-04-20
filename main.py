import re
import time
import requests
from bs4 import BeautifulSoup
import xml.etree.ElementTree as ET
import logging


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

file_handler = logging.FileHandler("blast.log")
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)


class Alignment:
    """
    A class that stores parsed BLAST alignment results, including
    full sequences and metadata.

    Parameters
    ----------
    subj_id : str
        Subject ID of the matched sequence.
    subj_name : str
        Subject sequence description.
    subj_len : int
        Length of the subject sequence.
    score_bits : float
        Score in bits from BLAST alignment.
    score : int
        Raw score from BLAST alignment.
    e_value : str
        E-value (expect value) from BLAST alignment.
    identities : int
        Number of identical matches.
    match_len : int
        Number of positions with aligned characters.
    align_len : int
        Total alignment length.
    query_align_chunks : list of tuple
        Query alignment segments as (start, sequence, end).
    sbjct_align_chunks : list of tuple
        Subject alignment segments as (start, sequence, end).
    """

    def __init__(
        self,
        subj_id,
        subj_name,
        subj_len,
        score_bits,
        score,
        e_value,
        identities,
        match_len,
        align_len,
        query_align_chunks,
        sbjct_align_chunks,
    ):
        self.subj_id = subj_id
        self.subj_name = subj_name
        self.subj_len = subj_len
        self.score = score
        self.score_bits = score_bits
        self.e_value = e_value
        self.match_len = match_len
        self.identities = identities
        self.align_len = align_len
        self.query_align_chunks = query_align_chunks
        self.sbjct_align_chunks = sbjct_align_chunks

        # полные строки выравнивания:
        self.query_align = "".join(seq for _, seq, _ in query_align_chunks)
        self.sbjct_align = "".join(seq for _, seq, _ in sbjct_align_chunks)

        self.subj_range = (sbjct_align_chunks[0][0], sbjct_align_chunks[-1][2])

    def __repr__(self):
        return (
            f"{self.subj_id=}, {self.subj_name=}, {self.subj_len=}, {self.subj_range=}, "
            f"{self.score=}, {self.e_value=}, {self.identities=}, {self.align_len=}, "
            f"{self.query_align=}, {self.sbjct_align=}"
        )


def run_blast(sequence, programm="tblastn", database="nt",
              taxon=None, wait=True, **params):
    """
    Submits a BLAST job to NCBI and retrieves parsed alignment results
    as Alignment objects.

    Parameters
    ----------
    sequence : str
        Protein or nucleotide sequence to be used as the query.
    programm : str
        Type of BLAST program to use (e.g., 'tblastn', 'blastn').
    database : str
        Database name or WGS prefixes to search against (e.g., 'wgs', 'nt',
        'WGS_VDB://XXXX01 WGS_VDB://XXXX02').
    taxon : str, optional
        Taxonomic restriction query string (e.g., species name or
        NCBI taxonomy ID).
    **params : dict
        Additional optional BLAST parameters.

    Returns
    -------
    alignments : list of Alignment
        Parsed results containing sequence alignments.
    """
    logger.info(f"Submitting BLAST: program={programm}, db={database}, taxon={taxon}")
    if database == "wgs" and taxon:
        response = requests.post(
            "https://www.ncbi.nlm.nih.gov/Traces/wgs/index.cgi?",
            data={'q': f'&wt=xml&q=text%3A*{taxon}*%20AND%20project_s%3Awgs'})
        prefixes = extract_prefix_organism_pairs(response.text)
        logger.info(f"prefixes: {prefixes}")
        logger.info(f"database: {database}")

        prefix_list = [p for p, _ in prefixes]
        prefixes = filter_valid_wgs_ids(prefix_list)  # checks for database validity
        for prefix, organism in prefixes.items():
            logger.info(f"{prefix}: {organism}")
        database = " ".join([f"WGS_VDB://{p}" for p in prefixes])
        logger.info(f"database: {database}")
        taxon = None

    data = {
        "CMD": "Put",
        "PROGRAM": programm,
        "DATABASE": database,
        "QUERY": sequence,
        "ENTREZ_QUERY": taxon,
    }

    # redefine values from kwargs params
    data.update(params)
    for key, value in params.items():
        logger.info(f"Extra param: {key} = {value}")

    # request
    response = requests.post("https://blast.ncbi.nlm.nih.gov/Blast.cgi", data=data)
    text = response.text
    logger.info(f"NCBI submission status: {response.status_code}")

    rid, rtoe = None, None
    for line in text.splitlines():
        if "RID =" in line:
            rid = line.split("=")[1].strip()
        if "RTOE =" in line:
            rtoe_str = line.split("=")[1].strip()
            try:
                rtoe = int(rtoe_str)
            except ValueError:
                rtoe = 20

    if rid:
        logger.info(f"RID: {rid} | Estimated wait: {rtoe} sec")
        logger.info(f"BLAST result link: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID={rid}")
    else:
        logger.error("RID not found in response")
        raise Exception("RID not found in response. Full response:\n" + text)

    if wait:
        return wait_for_blast_results(rid)
    else:
        return rid


def wait_for_blast_results(rid, rtoe=10, poll_interval=5, verbose=True):
    """
    Waits for a BLAST job to complete, fetches the result in text format,
    and parses the alignments.

    Parameters
    ----------
    rid : str
        Request ID (RID) returned from a previous BLAST submission.
    rtoe : int, optional
        Recommended time of execution (in seconds). Used as fallback sleep time.
    poll_interval : int, optional
        Interval in seconds between status checks while waiting.
    verbose : bool, optional
        Whether to logger.info status updates.

    Returns
    -------
    alignments : list of Alignment
        Parsed alignment results from the BLAST output.

    Raises
    ------
    Exception
        If the job fails, expires, or if result cannot be retrieved.
    """
    while True:
        response = requests.get(
            "https://blast.ncbi.nlm.nih.gov/Blast.cgi",
            params={"CMD": "Get", "RID": rid},
        )
        html = response.text

        if "There was a problem with the search" in html:
            logger.error(f"NCBI returned an error during check. RID: {rid}")
            raise Exception(f"NCBI returned an error during check. RID: {rid}")

        if "Status=WAITING" in html:
            if verbose:
                logger.info("Waiting for BLAST job to complete...")
            time.sleep(poll_interval)
        elif "Status=FAILED" in html:
            logger.error(f"BLAST search failed. RID: {rid}")
            raise Exception(f"BLAST search failed. RID: {rid}")
        elif "Status=UNKNOWN" in html:
            logger.error(f"BLAST RID expired or unknown. RID: {rid}")
            raise Exception(f"BLAST RID expired or unknown. RID: {rid}")
        elif "Status=READY" in html:
            if "dscTable" in html or "Sequences producing significant alignments" in html:
                break
            else:
                if verbose:
                    logger.info(
                        f"Warning: No significant hits found, but continuing to fetch result. RID: {rid}"
                    )
                break
        else:
            time.sleep(rtoe)

    result = requests.get(
        "https://blast.ncbi.nlm.nih.gov/Blast.cgi",
        params={
            "CMD": "Get",
            "FORMAT_OBJECT": "Alignment",
            "FORMAT_TYPE": "Text",
            "RID": rid,
            "DESCRIPTIONS": 100,
            "ALIGNMENTS": 100,
        },
    )

    if "There was a problem with the search" in result.text:
        logger.error(f"NCBI returned an error in final fetch. RID: {rid}")
        raise Exception(f"NCBI returned an error in final fetch. RID: {rid}")

    return parse_blast_text_output(result.text)


def parse_blast_text_output(text):
    """
    Parses BLAST text output and extracts all alignment blocks into Alignment objects.

    Parameters
    ----------
    text : str
        Raw BLAST output in text format.

    Returns
    -------
    alignments : list of Alignment
        List of alignment records parsed from the output.
    """
    alignments = []

    blocks = text.split("\n\n>")
    first_block = blocks[0]
    if first_block.startswith("<p>"):
        blocks[0] = first_block.split("ALIGNMENTS")[1].strip()
    else:
        blocks[0] = blocks[0]

    for block in blocks:
        lines = block.strip().split("\n")
        subj_line = lines[0]
        subj_id = subj_line.split()[0].lstrip(">")  # remove '>'
        subj_name = " ".join(subj_line.split()[1:])

        i = 1
        while i < len(lines):
            subj_len = score = e_value = identities = align_len = None
            query_chunks, sbjct_chunks = [], []

            while i < len(lines):
                line = lines[i]
                if line.startswith("Length="):
                    subj_len = int(line.split("=")[1])
                elif line.startswith(" Score ="):
                    score_match = re.search(
                        r"Score\s=\s([\d\.]+)\sbits\s\((\d+)\)", line
                    )
                    evalue_match = re.search(
                        r"Expect(?:\(\d+\))? = ([\deE\.\-]+)", line
                    )
                    score_bits = score_match.group(1) if score_match else None  # bits
                    score = score_match.group(2) if score_match else None
                    e_value = evalue_match.group(1) if evalue_match else None
                elif "Identities" in line:
                    ident_m = re.search(
                        r"Identities\s=\s(\d+)\/(\d+)\s\((\d+)\%\)", line
                    )
                    match_len = int(ident_m.group(1))
                    align_len = int(ident_m.group(2))
                    identities = int(ident_m.group(3))
                elif line.startswith("Query "):
                    query_parts = line.split()
                    q_start, q_seq, q_end = (
                        int(query_parts[1]),
                        query_parts[2],
                        int(query_parts[3]),
                    )
                    query_chunks.append((q_start, q_seq, q_end))

                    i += 2
                    if i < len(lines):
                        sbjct_parts = lines[i].split()
                        s_start, s_seq, s_end = (
                            int(sbjct_parts[1]),
                            sbjct_parts[2],
                            int(sbjct_parts[3]),
                        )
                        sbjct_chunks.append((s_start, s_seq, s_end))
                        # subj_range = (s_start, s_end)
                elif line.startswith(">") or line.startswith("Sequence ID:"):
                    break  # start of a new subject
                i += 1

            alignments.append(
                Alignment(
                    subj_id=subj_id,
                    subj_name=subj_name,
                    subj_len=subj_len,
                    score_bits=score_bits,
                    score=score,
                    e_value=e_value,
                    identities=identities,
                    align_len=align_len,
                    match_len=match_len,
                    query_align_chunks=query_chunks,
                    sbjct_align_chunks=sbjct_chunks,
                )
            )

            while i < len(lines) and not lines[i].startswith(" Score ="):
                i += 1
    logger.info(f"Parsed {len(alignments)} alignments from BLAST output")
    return alignments


def extract_prefix_organism_pairs(xml_text):
    root = ET.fromstring(xml_text)
    results = []

    # logger.info numFound
    result_element = root.find(".//result[@name='response']")
    if result_element is not None:
        num_found = result_element.attrib.get("numFound")
        logger.info(f"[WGS index] Found {num_found} matching entries.")

    for doc in root.findall(".//doc"):
        prefix = None
        organism = None
        for child in doc:
            if child.tag == "str":
                if child.attrib.get("name") == "prefix_s":
                    prefix = child.text
                elif child.attrib.get("name") == "organism_an":
                    organism = child.text
        if prefix and organism:
            results.append((prefix, organism))
    return results


def filter_valid_wgs_ids(prefixes):
    """
    Verifies WGS prefix validity through getDBInfo.cgi.

    Parameters
    ----------
    prefixes : list of str
        Prefix list (e.g., ['ACOL01', 'AEYK01'])

    Returns
    -------
    dict
        Dict {prefix: organism}, only for valid prefixes.
    """
    logger.info("Running filter_valid_wgs_ids")
    db_string = ",".join(f"WGS_VDB://{p}" for p in prefixes)

    headers = {
        "User-Agent": (
            "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
            "AppleWebKit/537.36 (KHTML, like Gecko) "
            "Chrome/134.0.0.0 Safari/537.36"),
        "Referer": "https://blast.ncbi.nlm.nih.gov/Blast.cgi",
        "Origin": "https://blast.ncbi.nlm.nih.gov",
        "Accept": "*/*",
    }

    params = {"DATABASE": db_string, "CMD": "getDBOrg"}

    response = requests.get(
        "https://blast.ncbi.nlm.nih.gov/getDBInfo.cgi", headers=headers, params=params
    )
    if response.status_code != 200:
        logger.error(f"Request failed with status code {response.status_code}")
        raise Exception(f"Request failed with status code {response.status_code}")

    soup = BeautifulSoup(response.text, "html.parser")
    table = soup.find("table", {"id": "dbSpecies"})
    if not table:
        logger.error(f"No species table found in response")
        raise Exception("No species table found in response.")

    valid = {}
    for row in table.find_all("tr")[1:]:
        cols = row.find_all("td")
        if len(cols) >= 2:
            db = cols[0].text.strip().replace("WGS_VDB://", "")
            organism = cols[1].text.strip()
            valid[db] = organism

    logger.info(f"Validated WGS prefixes: {list(valid.keys())}")
    return valid
