#!/usr/bin/env python3
"""Query ADS using an ephemeral, hidden prompt; save only public responses."""
import argparse
import getpass
import json
from pathlib import Path
import requests
from adastra.archive import sha256_file, utc_now, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--queries', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() and any(args.output.iterdir()):
        raise ValueError('Refusing to overwrite an existing literature search')
    queries = json.loads(args.queries.read_text())
    token = getpass.getpass('ADS token (hidden, used only in this process): ')
    if not token.strip():
        raise ValueError('Missing ADS token')
    declaration = {'created_utc': utc_now(), 'queries_path': str(args.queries),
                   'queries_sha256': sha256_file(args.queries),
                   'source_sha256': sha256_file(Path(__file__)),
                   'endpoint': 'https://api.adsabs.harvard.edu/v1/search/query',
                   'scope': 'Scoped identity and topic queries; a nonmatch is not novelty clearance',
                   'credentials_saved': False, 'retry_policy': 'Single attempt per query'}
    write_json(args.output / 'declaration.json', declaration)
    summaries = []
    for item in queries:
        params = {'q': item['q'], 'rows': 100, 'start': 0,
                  'fl': 'bibcode,title,author,year,doi,abstract,identifier', 'sort': 'date desc'}
        row = {'name': item['name'], 'query': params, 'retrieved_utc': utc_now()}
        try:
            response = requests.get(declaration['endpoint'], params=params,
                                    headers={'Authorization': 'Bearer ' + token.strip()}, timeout=60)
            row['http_status'] = response.status_code
            # Error bodies and request headers are deliberately not persisted.
            response.raise_for_status()
            payload = response.json()
            found = payload['response']['numFound']
            docs = payload['response']['docs']
            path = args.output / (item['name'] + '.json')
            write_json(path, payload)
            row.update(status='complete' if found <= len(docs) else 'truncated',
                       num_found=found, returned=len(docs), raw_path=str(path), sha256=sha256_file(path))
        except requests.RequestException as error:
            row.update(status='failed', error_type=type(error).__name__)
        except (ValueError, KeyError) as error:
            row.update(status='failed', error_type=type(error).__name__)
        summaries.append(row)
        write_json(args.output / 'summary.json', {'declaration': declaration, 'queries': summaries,
                                                  'completed_utc': utc_now()})
        print(f"{row['name']}: {row['status']} HTTP {row.get('http_status')} count {row.get('num_found')}", flush=True)
        if row.get('http_status') in (401, 403):
            break


if __name__ == '__main__':
    main()
