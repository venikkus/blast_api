import re
import time
import requests
from bs4 import BeautifulSoup
import xml.etree.ElementTree as ET
import logging

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


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
              taxon=None, taxids=None, exclude_taxids=None, wait=True, **params):
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
    taxids : iterable of int or str, optional
        NCBI taxonomy IDs whose WGS projects should be searched. This uses
        NCBI's taxid2wgs endpoint and accepts higher-level taxonomic groups.
    exclude_taxids : iterable of int or str, optional
        Taxonomy IDs to exclude when resolving WGS projects.
    **params : dict
        Additional optional BLAST parameters.

    Returns
    -------
    alignments : list of Alignment
        Parsed results containing sequence alignments.
    """
    if database == "wgs" and (taxids is not None or taxon is not None):
        if taxids is not None and taxon is not None:
            raise ValueError("Pass either taxon or taxids, not both")
        if taxids is not None:
            include_taxids = taxids
        elif str(taxon).isdigit():
            include_taxids = [taxon]
        else:
            include_taxids = [resolve_taxon_id(taxon)]
        projects = get_wgs_projects(include_taxids, exclude_taxids)
        if not projects:
            raise ValueError("NCBI returned no WGS projects for the selected taxids")
        database = " ".join(projects)
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
        print(f"{key}: {str(value)}")

    # request
    response = requests.post("https://blast.ncbi.nlm.nih.gov/Blast.cgi", data=data)
    text = response.text
    print(f"NCBI initial response: {response.status_code}")

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
        print(f"RID: {rid} | Estimated wait: {rtoe} sec")
        print(f"BLAST result link: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID={rid}")
    else:
        raise Exception("RID not found in response. Full response:\n" + text)

    if wait:
        return wait_for_blast_results(rid)
    else:
        return rid


def wait_for_blast_results(rid, rtoe=10, poll_interval=5, verbose=True,
                           alignment_limit=100):
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
        Whether to print status updates.

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
            raise Exception(f"NCBI returned an error during check. RID: {rid}")

        if "Status=WAITING" in html:
            if verbose:
                print("Waiting for BLAST job to complete...")
            time.sleep(poll_interval)
        elif "Status=FAILED" in html:
            raise Exception(f"BLAST search failed. RID: {rid}")
        elif "Status=UNKNOWN" in html:
            raise Exception(f"BLAST RID expired or unknown. RID: {rid}")
        elif "Status=READY" in html:
            if "dscTable" in html or "Sequences producing significant alignments" in html:
                break
            else:
                if verbose:
                    print(
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
            "DESCRIPTIONS": alignment_limit,
            "ALIGNMENTS": alignment_limit,
        },
    )

    if "There was a problem with the search" in result.text:
        raise Exception(f"NCBI returned an error in final fetch. RID: {rid}")

    if "No hits found" in result.text or "No significant similarity found" in result.text:
        return []
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

    return alignments


def get_wgs_projects(taxids, exclude_taxids=None, timeout=60):
    """Resolve taxonomy IDs to NCBI WGS_VDB project identifiers.

    Uses the endpoint behind NCBI's official ``taxid2wgs.pl`` utility.
    ``taxids`` may contain one or more species or higher-level taxonomy IDs.
    """
    if isinstance(taxids, (str, int)):
        taxids = [taxids]
    include_ids = [str(taxid).strip() for taxid in taxids]
    include_ids = [taxid for taxid in include_ids if taxid]
    if not include_ids or any(not taxid.isdigit() for taxid in include_ids):
        raise ValueError("taxids must contain one or more numeric NCBI TaxIDs")

    params = {"INCLUDE_TAXIDS": ",".join(include_ids)}
    if exclude_taxids is not None:
        if isinstance(exclude_taxids, (str, int)):
            exclude_taxids = [exclude_taxids]
        exclude_ids = [str(taxid).strip() for taxid in exclude_taxids]
        exclude_ids = [taxid for taxid in exclude_ids if taxid]
        if any(not taxid.isdigit() for taxid in exclude_ids):
            raise ValueError("exclude_taxids must contain numeric NCBI TaxIDs")
        if exclude_ids:
            params["EXCLUDE_TAXIDS"] = ",".join(exclude_ids)

    response = requests.post(
        "https://www.ncbi.nlm.nih.gov/blast/BDB2EZ/taxid2wgs.cgi",
        params=params,
        timeout=timeout,
    )
    response.raise_for_status()
    projects = []
    for token in response.text.split():
        if token.startswith("WGS_VDB://"):
            project = token[len("WGS_VDB://"):]
            if re.fullmatch(r"[A-Z]{4,6}[0-9]{2}", project):
                projects.append(token)
    return list(dict.fromkeys(projects))


def resolve_taxon_id(name, timeout=30):
    """Resolve one exact NCBI Taxonomy name or synonym to a TaxID.

    Ambiguous names raise ``ValueError`` with candidate taxa; callers can then
    pass the intended identifier explicitly via ``taxids``.
    """
    name = str(name).strip()
    if not name:
        raise ValueError("Taxon name cannot be empty")
    if name.isdigit():
        return name

    response = requests.get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
        params={"db": "taxonomy", "term": f'"{name}"[name]',
                "retmode": "xml", "retmax": 20},
        timeout=timeout,
    )
    response.raise_for_status()
    root = ET.fromstring(response.content)
    taxids = [element.text for element in root.findall("./IdList/Id") if element.text]

    if not taxids:
        raise LookupError(f"No NCBI Taxonomy match found for {name!r}")
    if len(taxids) > 1:
        summary_response = requests.get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
            params={"db": "taxonomy", "id": ",".join(taxids), "retmode": "json"},
            timeout=timeout,
        )
        summary_response.raise_for_status()
        summaries = summary_response.json().get("result", {})
        candidates = []
        for taxid in taxids:
            record = summaries.get(taxid, {})
            taxname = record.get("taxname", "unknown name")
            rank = record.get("rank", "unknown rank")
            candidates.append(f"{taxname} ({rank}, TaxID {taxid})")
        raise ValueError(
            f"Taxon name {name!r} is ambiguous: " + "; ".join(candidates)
            + ". Pass the intended NCBI TaxID via taxids."
        )
    return taxids[0]


def extract_prefix_organism_pairs(xml_text):
    root = ET.fromstring(xml_text)
    results = []

    # print numFound
    result_element = root.find(".//result[@name='response']")
    if result_element is not None:
        num_found = result_element.attrib.get("numFound")
        print(f"[WGS index] Found {num_found} matching entries.")

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


def filter_valid_wgs_ids(prefixes, batch_size=50):
    """
    Verifies WGS prefix validity through getDBInfo.cgi.

    Parameters
    ----------
    prefixes : list of str
        Prefix list (e.g., ['ACOL01', 'AEYK01'])
    batch_size : int
        How many prefixes to check at once to avoid 414 URI Too Long

    Returns
    -------
    dict
        Dict {prefix: organism}, only for valid prefixes.
    """
    db_string = ",".join(f"WGS_VDB://{p}" for p in prefixes)
    all_valid = {}

    headers = {
        "User-Agent": (
            "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
            "AppleWebKit/537.36 (KHTML, like Gecko) "
            "Chrome/134.0.0.0 Safari/537.36"),
        "Referer": "https://blast.ncbi.nlm.nih.gov/Blast.cgi",
        "Origin": "https://blast.ncbi.nlm.nih.gov",
        "Accept": "*/*",
    }

    for i in range(0, len(prefixes), batch_size):
        batch = prefixes[i:i + batch_size]
        db_string = ",".join(f"WGS_VDB://{p}" for p in batch)
        params = {"DATABASE": db_string, "CMD": "getDBOrg"}

        logger.info(f"Filtering batch {i//batch_size + 1}: {len(batch)} prefixes")

        response = requests.get(
            "https://blast.ncbi.nlm.nih.gov/getDBInfo.cgi",
            headers=headers,
            params=params
        )
        if response.status_code != 200:
            logger.error(f"Request failed with status code {response.status_code}")
            raise Exception(f"Request failed with status code {response.status_code}")

        soup = BeautifulSoup(response.text, "html.parser")
        table = soup.find("table", {"id": "dbSpecies"})
        if not table:
            logger.warning("No species table found in response for this batch")
            continue

        for row in table.find_all("tr")[1:]:
            cols = row.find_all("td")
            if len(cols) >= 2:
                db = cols[0].text.strip().replace("WGS_VDB://", "")
                organism = cols[1].text.strip()
                all_valid[db] = organism

    logger.info(f"Total validated prefixes: {len(all_valid)}")
    return all_valid
