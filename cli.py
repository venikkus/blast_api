import argparse
import ast
from main import run_blast, extract_prefix_organism_pairs, filter_valid_wgs_ids
from pathlib import Path
import logging
import requests

logging.basicConfig(
    filename='blast_cli.log',
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)


def parse_input(input_arg):
    path = Path(input_arg)
    if path.exists():
        if path.is_dir():
            for filename in path.iterdir():
                if filename.suffix in ('.fasta', '.fa', '.txt'):
                    yield from parse_fasta(filename)
        elif path.is_file():
            yield from parse_fasta(path)
    else:
        yield ("query_from_cli", input_arg)


def parse_fasta(file_path):
    sequences = []
    with open(file_path, 'r', encoding='utf-8') as file:
        content = file.read()
        entries = content.strip().split('>')
        for entry in entries:
            if entry:
                lines = entry.strip().split('\n')
                header = lines[0]
                seq = ''.join(lines[1:])
                sequences.append((header, seq))
    return sequences


def cli_run_blast(args):
    extra_params = {}
    if args.extra:
        for i in range(0, len(args.extra), 2):
            key = args.extra[i].lstrip('-')
            value = args.extra[i + 1]
            extra_params[key] = value

    for header, sequence in parse_input(args.input):
        print(f"Processing {header}")
        alignments = run_blast(
            sequence=sequence,
            programm=args.program,
            database=args.database,
            taxon=args.taxon,
            **extra_params
        )
        for aln in alignments:
            print(aln)


def cli_extract_prefixes(args):
    response = requests.post('https://www.ncbi.nlm.nih.gov/Traces/wgs/index.cgi?',
                             data={'q': f'&wt=xml&rows={rows}&sort=prefix_s%20asc&q=text%3A*{args.taxon}*%20AND%20project_s%3Awgs&stats=true'})
    prefixes = extract_prefix_organism_pairs(response.text)
    for prefix, organism in prefixes:
        print(f"{prefix}: {organism}")


def cli_validate_prefixes(args):
    if args.prefixes:
        try:
            prefixes = ast.literal_eval(args.prefixes)
            if not isinstance(prefixes, list):
                raise ValueError
        except Exception:
            raise ValueError("Argument --prefixes must be a string representation of a list, e.g., \"['ACOL01']\"")
    else:
        import requests
        response = requests.post(
            "https://www.ncbi.nlm.nih.gov/Traces/wgs/index.cgi?",
            data={'q': f'&wt=xml&q=text%3A*{args.taxon}*%20AND%20project_s%3Awgs'}
        )
        pairs = extract_prefix_organism_pairs(response.text)
        prefixes = [p for p, _ in pairs]

    valid = filter_valid_wgs_ids(prefixes)

    print(f"Validation summary: {len(prefixes)} → {len(valid)}")
    for prefix, organism in valid.items():
        print(f"{prefix}: {organism}")


def cli():
    parser = argparse.ArgumentParser(prog='blast-cli', description='Unified BLAST CLI tool')
    subparsers = parser.add_subparsers(dest='command')

    # blast
    parser_blast = subparsers.add_parser('blast', help='Run BLAST for all sequences in directory')
    parser_blast.add_argument('-i', '--input', required=True,
                              help='FASTA file, folder with FASTA files, or raw sequence string')
    parser_blast.add_argument('-p', '--program', default='blastn', help='BLAST program')
    parser_blast.add_argument('-db', '--database', default='nt', help='Database (nt, wgs, etc)')
    parser_blast.add_argument('-t', '--taxon', help='Taxon name')
    parser_blast.add_argument('--extra', nargs=argparse.REMAINDER, help='Extra BLAST params')
    parser_blast.set_defaults(func=cli_run_blast)

    # extract
    parser_extract = subparsers.add_parser('extract', help='Extract WGS prefixes from NCBI by taxon')
    parser_extract.add_argument('-t', '--taxon', required=True, help='Taxon')
    parser_extract.set_defaults(func=cli_extract_prefixes)

    # validate
    parser_validate = subparsers.add_parser('validate', help='Validate WGS prefixes via NCBI')
    group = parser_validate.add_mutually_exclusive_group(required=True)
    group.add_argument('-p', '--prefixes', help="List of prefixes, e.g., \"['ACOL01']\"")
    group.add_argument('-t', '--taxon', help='Taxon name')
    parser_validate.set_defaults(func=cli_validate_prefixes)

    args = parser.parse_args()
    if hasattr(args, 'func'):
        args.func(args)
    else:
        parser.print_help()


if __name__ == '__main__':
    cli()
