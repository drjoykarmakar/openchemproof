#!/usr/bin/env python3
"""Recreate a hash-pinned 100-record ORD evidence subset.

Dependencies: ord-schema==0.9.6, pyarrow==26.0.1, requests==2.34.2.
Run with a separate environment; this does not modify OpenChemProof dependencies.
"""
from __future__ import annotations
import argparse
import collections
import hashlib
import json
from pathlib import Path
import tempfile

import pyarrow.parquet as pq
import requests
from google.protobuf.json_format import MessageToDict
from ord_schema.proto import reaction_pb2

ROOT = Path(__file__).resolve().parent
GITHUB_COMMIT = '398d7349d4218eccf8782a350a13588c7c7fb88f'
HF_COMMIT = '69a17d69acd264c1ccde52c1a86ef7689c43a7c2'
SOURCES = [
    ('ord_dataset-99c23cf435dc42f1af884053bc8b11c7', '59cd8c35e021be300944f8e556b8fd7ec30a6bc32ff36a54be133aad8e44c393', 214444),
    ('ord_dataset-00005539a1e04c809a9a78647bea649c', '5326c1f79ef6a71e6b2a91cbd07cd8521049427df79307686c7c4f44eb16b53a', 285104),
]


def canonical(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(',', ':'))


def sha(data):
    return hashlib.sha256(data).hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, sort_keys=True, indent=2) + '\n')


def write_jsonl(path, values):
    path.write_text(''.join(canonical(v) + '\n' for v in values))


def yields(product):
    # Protobuf defaults are deliberately not used: absent value is not zero.
    return [m['percentage']['value'] for m in product.get('measurements', [])
            if m.get('type') == 'YIELD' and 'value' in m.get('percentage', {})]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache-dir', type=Path, default=Path(tempfile.gettempdir()) / 'openchemproof-ord-scan')
    parser.add_argument('--offline', action='store_true', help='Only use already downloaded, hash-verified source files')
    args = parser.parse_args()
    args.cache_dir.mkdir(parents=True, exist_ok=True)
    selection = json.loads((ROOT / 'selection.json').read_text())
    assert len(selection) == 100 and len({s['reaction_id'] for s in selection}) == 100
    wanted = {s['reaction_id'] for s in selection}
    originals, protobuf_hashes, source_lookup, manifests = {}, {}, {}, []
    source_counts = {}
    for dataset_id, expected_hash, expected_size in SOURCES:
        filename = dataset_id + '.parquet'
        relative = f'data/{dataset_id[12:14]}/{filename}'
        download_url = f'https://huggingface.co/datasets/open-reaction-database/ord-data/resolve/{HF_COMMIT}/{relative}'
        cache_path = args.cache_dir / filename
        if not cache_path.exists():
            if args.offline:
                raise FileNotFoundError(cache_path)
            response = requests.get(download_url, timeout=120)
            response.raise_for_status()
            payload = response.content
            if sha(payload) != expected_hash or len(payload) != expected_size:
                raise ValueError('Source download hash/size mismatch: ' + dataset_id)
            cache_path.write_bytes(payload)
        payload = cache_path.read_bytes()
        if sha(payload) != expected_hash or len(payload) != expected_size:
            raise ValueError('Cached source hash/size mismatch: ' + dataset_id)
        table = pq.read_table(cache_path)
        meta = {k.decode(): v.decode() for k, v in table.schema.metadata.items() if k.startswith(b'ord.')}
        assert meta['ord.dataset_id'] == dataset_id
        manifest = {
            'dataset_id': dataset_id, 'dataset_name': meta['ord.name'],
            'dataset_description': meta.get('ord.description'), 'source_reaction_count': table.num_rows,
            'github_commit': GITHUB_COMMIT, 'huggingface_commit': HF_COMMIT,
            'source_url': f'https://github.com/open-reaction-database/ord-data/blob/{GITHUB_COMMIT}/{relative}',
            'download_url': download_url, 'sha256': expected_hash, 'bytes': expected_size,
            'license': 'CC-BY-SA-4.0',
            'license_url': f'https://github.com/open-reaction-database/ord-data/blob/{GITHUB_COMMIT}/LICENSE',
        }
        manifests.append(manifest)
        counters = collections.Counter()
        for row in table.to_pylist():
            proto = reaction_pb2.Reaction.FromString(row['reaction'])
            record = MessageToDict(proto, preserving_proto_field_name=True)
            assert row['reaction_id'] == record['reaction_id']
            all_yields = [v for o in record.get('outcomes', []) for p in o.get('products', []) for v in yields(p)]
            counters['reactions_with_explicit_zero_yield'] += 0 in all_yields
            counters['reactions_with_explicit_positive_yield'] += any(v > 0 for v in all_yields)
            if record['reaction_id'] in wanted:
                originals[record['reaction_id']] = record
                protobuf_hashes[record['reaction_id']] = sha(row['reaction'])
                source_lookup[record['reaction_id']] = manifest
        source_counts[dataset_id] = dict(counters)
    assert wanted == set(originals), 'Selected record missing from pinned source'
    normalized, raw = [], []
    for chosen in selection:
        rid = chosen['reaction_id']
        record, source = originals[rid], source_lookup[rid]
        assert source['dataset_id'] == chosen['dataset_id']
        oi, pi = chosen['outcome_index'], chosen['product_index']
        outcome = record['outcomes'][oi]
        product = outcome['products'][pi]
        values = yields(product)
        identity_keys = [m['analysis_key'] for m in product.get('measurements', []) if m.get('type') == 'IDENTITY' and m.get('analysis_key')]
        identity_analyses = {key: outcome['analyses'][key] for key in identity_keys}
        is_success = chosen['cohort'] == 'source_characterized_success'
        if is_success:
            assert any(50 <= v <= 100 for v in values)
            types = {a.get('type') for a in identity_analyses.values() if a.get('data')}
            assert {'NMR_1H', 'HRMS'} <= types
        else:
            assert values == [0.0], (rid, values)
            for excerpt in chosen['source_failure_excerpts']:
                assert excerpt in record['notes']['procedure_details']
        raw.append(record)
        normalized.append({
            'record_id': rid, 'schema_version': 'openchemproof.ord-evidence.v1',
            'cohort': chosen['cohort'], 'source_dataset_id': source['dataset_id'],
            'source_dataset_name': source['dataset_name'], 'source_url': source['source_url'],
            'source_download_url': source['download_url'], 'source_dataset_sha256': source['sha256'],
            'source_reaction_protobuf_sha256': protobuf_hashes[rid],
            'source_reaction_canonical_json_sha256': sha(canonical(record).encode()),
            'source_record_locator': {'parquet_column': 'reaction', 'reaction_id': rid, 'outcome_index': oi, 'product_index': pi},
            'license': 'CC-BY-SA-4.0', 'attribution': source['dataset_name'] + '; Open Reaction Database; see source provenance and sources_ord.md',
            'publication_doi': record.get('provenance', {}).get('doi'),
            'publication_url': record.get('provenance', {}).get('publication_url'),
            'inputs': record.get('inputs', {}), 'conditions': record.get('conditions', {}),
            'procedure_details': record.get('notes', {}).get('procedure_details'),
            'product_identifiers': product.get('identifiers', []),
            'recorded_yield_percentages': values,
            'yield_measurements': [m for m in product.get('measurements', []) if m.get('type') == 'YIELD'],
            'yield_analyses': {m['analysis_key']: outcome.get('analyses', {}).get(m['analysis_key']) for m in product.get('measurements', []) if m.get('type') == 'YIELD' and m.get('analysis_key')},
            'identity_analyses': identity_analyses,
            'characterization_evidence': 'product-linked NMR peak text and observed HRMS m/z' if is_success else 'failure narrative and recorded zero yield',
            'source_failure_excerpts': chosen.get('source_failure_excerpts', []),
            'review_note': chosen['review_note'],
            'raw_spectral_traces_available': False,
            'independently_verified_spectrum_assignment': False,
            'engine_confirmed_result': None,
            'strict_negative_eligible': False,
            'strict_negative_missing': [] if is_success else ['validated detection limit', 'documented positive/negative controls', 'complete assay trace and target-specific calibration'],
        })
    (ROOT / 'raw').mkdir(exist_ok=True)
    write_jsonl(ROOT / 'raw' / 'selected_reactions.jsonl', raw)
    write_jsonl(ROOT / 'records.jsonl', normalized)
    summary = {
        'total_records': len(normalized), 'source_characterized_successes': sum(r['cohort'] == 'source_characterized_success' for r in normalized),
        'source_reported_zero_yield_with_failure_narrative': sum(r['cohort'] == 'source_reported_zero_yield' for r in normalized),
        'successes_yield_at_least_50_percent': sum(r['cohort'] == 'source_characterized_success' and max(r['recorded_yield_percentages']) >= 50 for r in normalized),
        'successes_with_product_linked_nmr_peak_text_and_hrms': 60,
        'complete_raw_spectra_cases': 0, 'independently_verified_spectral_match_cases': 0,
        'engine_confirmed_positive_cases': 0, 'strict_negative_eligible_cases': 0,
        'target_shortfalls': {'source_characterized_successes_of_60': 0, 'source_reported_zero_yield_of_40': 0, 'raw_spectra_verified_successes_of_60': 60, 'strict_controlled_negatives_of_40': 40},
        'positive_yield_range_percent': [min(min(r['recorded_yield_percentages']) for r in normalized if r['cohort'] == 'source_characterized_success'), max(max(r['recorded_yield_percentages']) for r in normalized if r['cohort'] == 'source_characterized_success')],
        'success_publication_dois': sorted({r['publication_doi'] for r in normalized if r['cohort'] == 'source_characterized_success'}),
        'pinned_source_dataset_counts': source_counts,
        'selection_bias': 'Curated provenance demonstration, not a random or balanced general-chemistry benchmark. Success borylations and failed Buchwald-Hartwig reactions differ in source and chemistry.',
    }
    write_json(ROOT / 'summary.json', summary)
    write_json(ROOT / 'manifest.json', {
        'acquired_date_utc': '2026-11-16', 'source_license': 'CC-BY-SA-4.0',
        'sources': manifests, 'selection_file': 'selection.json', 'selection_method': 'First eligible success records round-robin by publication DOI; zero-yield records selected after narrative review. Selection IDs frozen.',
        'dependencies': {'ord-schema': '0.9.6', 'pyarrow': '26.0.1', 'requests': '2.34.2'},
        'outputs': {str(p.relative_to(ROOT)): {'bytes': p.stat().st_size, 'sha256': sha(p.read_bytes())} for p in [ROOT/'raw'/'selected_reactions.jsonl', ROOT/'records.jsonl', ROOT/'selection.json', ROOT/'summary.json', ROOT/'LICENSE-CC-BY-SA-4.0.txt']},
    })
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
