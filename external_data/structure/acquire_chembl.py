#!/usr/bin/env python3
"""Acquire exact-assay stereoisomer pairs from the public ChEMBL API.

Data and derivative records: CC-BY-SA-3.0; this acquisition code follows
OpenChemProof's code license. Requires Python 3.11+ and RDKit.
"""
from __future__ import annotations
import argparse, concurrent.futures, hashlib, itertools, json, time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
from rdkit import Chem, rdBase

ROOT = Path(__file__).resolve().parent
BASE = 'https://www.ebi.ac.uk/chembl/api/data/'
TARGET = 'CHEMBL217'
KEEP = ['activity_id','assay_chembl_id','assay_description','assay_type','assay_variant_accession','assay_variant_mutation','activity_properties','canonical_smiles','data_validity_comment','data_validity_description','document_chembl_id','document_journal','document_year','molecule_chembl_id','molecule_pref_name','parent_molecule_chembl_id','pchembl_value','potential_duplicate','relation','src_id','standard_flag','standard_relation','standard_type','standard_units','standard_value','target_chembl_id','target_organism','target_pref_name','type','units','value','activity_comment']

def request(resource, params=None):
    url = BASE + resource + (('?' + urlencode(params)) if params else '')
    for attempt in range(4):
        try:
            with urlopen(url, timeout=180) as response:
                return json.load(response), url
        except Exception:
            if attempt == 3: raise
            time.sleep(2 ** attempt)

def acquire_page(offset):
    path = ROOT / 'raw' / f'activity_{offset:06d}.json'
    if path.exists(): return json.loads(path.read_text())
    params = {'target_chembl_id':TARGET,'standard_type__in':'Ki,IC50','standard_relation':'=','standard_units':'nM','order_by':'activity_id','limit':1000,'offset':offset}
    body, url = request('activity.json', params)
    result = {'retrieved_at':datetime.now(timezone.utc).isoformat(), 'url':url, 'page_meta':body['page_meta'], 'activities':[{key:a.get(key) for key in KEEP} for a in body['activities']]}
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2)+'\n')
    print(f'Fetched {offset}: {len(result["activities"])}', flush=True)
    return result

def structure(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None: return None
    chiral = Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)
    # Require explicitly resolved tetrahedral stereochemistry; unknown stereo
    # and racemate-versus-enantiomer comparisons are not stereo pairs here.
    if not chiral or any(label == '?' for _,label in chiral): return None
    iso = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    plain = Chem.Mol(mol); Chem.RemoveStereochemistry(plain)
    return {'canonical_isomeric_smiles':iso,'connectivity_smiles':Chem.MolToSmiles(plain,canonical=True),'rdkit_inchikey':Chem.MolToInchiKey(mol),'specified_tetrahedral_centers':len(chiral)}

def candidates(activities):
    groups = defaultdict(list); cache={}
    for a in activities:
        if a['data_validity_comment'] or a['potential_duplicate'] or a['standard_relation']!='=' or not a['standard_value'] or float(a['standard_value'])<=0: continue
        smi=a.get('canonical_smiles')
        if not smi: continue
        if smi not in cache: cache[smi]=structure(smi)
        s=cache[smi]
        if not s: continue
        k=(s['connectivity_smiles'],a['assay_chembl_id'],a['document_chembl_id'],a['standard_type'],a['standard_units'],a['target_chembl_id'])
        groups[k].append((a,s))
    result=[]; seen=set()
    for k,group in sorted(groups.items(),key=lambda x:min(a['activity_id'] for a,s in x[1])):
        by_iso={}
        for a,s in sorted(group,key=lambda x:x[0]['activity_id']):
            by_iso.setdefault(s['canonical_isomeric_smiles'],(a,s))
        for left,right in itertools.combinations(by_iso.values(),2):
            a,b=left[0],right[0]
            # Select distinct molecules only once across the complete fixture.
            if a['molecule_chembl_id'] in seen or b['molecule_chembl_id'] in seen: continue
            if any(a.get(field)!=b.get(field) for field in ['assay_description','assay_variant_accession','assay_variant_mutation','activity_properties','target_organism']): continue
            result.append((left,right)); seen.update([a['molecule_chembl_id'],b['molecule_chembl_id']])
    return result

def main():
    p=argparse.ArgumentParser();p.add_argument('--max-pages',type=int,default=10);p.add_argument('--pairs',type=int,default=50);p.add_argument('--offline',action='store_true');args=p.parse_args()
    if not args.offline:
        status,status_url=request('status.json')
        status_path=ROOT/'api_status.json'
        if status_path.exists():
            cached_status=json.loads(status_path.read_text())['status']
            if cached_status['chembl_db_version']!=status['chembl_db_version']:
                raise RuntimeError('ChEMBL release differs from cache; use --offline or a fresh cache directory.')
        else:
            status_path.write_text(json.dumps({'url':status_url,'retrieved_at':datetime.now(timezone.utc).isoformat(),'status':status},indent=2)+'\n')
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool: list(pool.map(acquire_page,range(0,args.max_pages*1000,1000)))
    pages=[json.loads(f.read_text()) for f in sorted((ROOT/'raw').glob('activity_*.json'))]
    activities=[a for page in pages for a in page['activities']]
    pairs=candidates(activities)[:args.pairs]
    print(f'Candidate pairs selected: {len(pairs)} from {len(activities)} activity records',flush=True)
    ids={'assay':sorted({x[0]['assay_chembl_id'] for pair in pairs for x in pair}), 'document':sorted({x[0]['document_chembl_id'] for pair in pairs for x in pair})}
    def load_metadata(kind):
        path=ROOT/'raw'/f'{kind}_metadata.json'
        if path.exists():
            cached=json.loads(path.read_text())
            available={item[kind+'_chembl_id'] for item in cached['response'][kind+'s']}
            if set(ids[kind])<=available: return cached
        if args.offline: raise RuntimeError(f'Missing metadata cache: {path}')
        body,url=request(kind+'.json',{kind+'_chembl_id__in':','.join(ids[kind]),'limit':1000})
        result={'url':url,'retrieved_at':datetime.now(timezone.utc).isoformat(),'response':body}
        path.write_text(json.dumps(result,indent=2)+'\n')
        return result
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        loaded=dict(zip(ids,pool.map(load_metadata,ids)))
    assays={a['assay_chembl_id']:a for a in loaded['assay']['response']['assays']}
    documents={d['document_chembl_id']:d for d in loaded['document']['response']['documents']}
    rows=[]; pair_rows=[]
    for n,(left,right) in enumerate(pairs,1):
        pair_id=f'chembl_stereo_{n:03d}'
        values=[float(x[0]['standard_value']) for x in (left,right)]
        ratio=max(values)/min(values)
        for side,(a,s) in zip('AB',(left,right)):
            rows.append({'record_id':pair_id+'_'+side,'pair_id':pair_id,'source':'ChEMBL','source_record_id':str(a['activity_id']),'source_url':BASE+f'activity/{a["activity_id"]}.json','molecule_url':'https://www.ebi.ac.uk/chembl/explore/compound/'+a['molecule_chembl_id'],'assay_url':BASE+f'assay/{a["assay_chembl_id"]}.json','document_url':BASE+f'document/{a["document_chembl_id"]}.json','source_license':'CC-BY-SA-3.0','case_type':'stereoisomer_same_assay','evidence_class':'assay_matched_numeric_activity','structure_identity':s,'measurement':a,'assay_metadata':assays[a['assay_chembl_id']],'document_metadata':documents[a['document_chembl_id']],'condition_completeness':'Source assay description and activity_properties retained verbatim; unreported conditions are unknown.','aggregation_status':'not_assessed','tautomer_status':'not_assessed','limitations':['Same ChEMBL assay identifier is used as the comparability criterion; complete experimental conditions and uncertainty are not necessarily recorded.','Numeric potency differences are observed values, not proof of a statistically significant or binary activity flip.']})
        pair_rows.append({'pair_id':pair_id,'record_ids':[pair_id+'_A',pair_id+'_B'],'molecule_chembl_ids':[left[0]['molecule_chembl_id'],right[0]['molecule_chembl_id']],'assay_chembl_id':left[0]['assay_chembl_id'],'endpoint':left[0]['standard_type'],'units':left[0]['standard_units'],'values':values,'observed_potency_ratio':ratio,'connectivity_equal':True,'isomeric_smiles_differ':True,'activity_flip_label':None,'comparison_scope':'same_assay_same_document_same_endpoint_same_units','uncertainty':'No significance test or binary flip threshold applied.'})
    for name,data in [('records.jsonl',rows),('pairs.jsonl',pair_rows)]:
        (ROOT/name).write_text(''.join(json.dumps(r,sort_keys=True)+'\n' for r in data))
    manifest={'dataset':'OpenChemProof ChEMBL stereoisomer evidence subset','created_at':datetime.now(timezone.utc).isoformat(),'source':'ChEMBL','source_homepage':'https://www.ebi.ac.uk/chembl/','license':'CC-BY-SA-3.0','license_url':'https://creativecommons.org/licenses/by-sa/3.0/','license_evidence':'https://chembl.gitbook.io/chembl-interface-documentation/frequently-asked-questions/general-questions','source_citation':'Zdrazil et al. The ChEMBL Database in 2023. Nucleic Acids Research. DOI:10.1093/nar/gkad1004','records':len(rows),'pairs':len(pairs),'unique_molecules':len({r['measurement']['molecule_chembl_id'] for r in rows}),'unique_assays':len({r['measurement']['assay_chembl_id'] for r in rows}),'screened_activities':len(activities),'target_id':TARGET,'rdkit_version':rdBase.rdkitVersion,'confirmed_aggregator_records':0,'tautomer_records':0,'api_status':json.loads((ROOT/'api_status.json').read_text()) if (ROOT/'api_status.json').exists() else None,'selection':'First qualifying non-overlapping fully assigned tetrahedral stereoisomer pairs in activity_id order, grouped by stereo-stripped exact canonical connectivity, same assay, document, endpoint, units, description, assay variant and recorded activity properties. Validity warnings and potential duplicate flags excluded. No activity ratio threshold.','raw_fields_excluded':['ligand_efficiency','commercial calculated molecule properties'],'files':{str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in [ROOT/'records.jsonl',ROOT/'pairs.jsonl',ROOT/'acquire_chembl.py',ROOT/'api_status.json',*sorted((ROOT/'raw').glob('*.json'))]}}
    (ROOT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({k:manifest[k] for k in ('records','pairs','unique_molecules','unique_assays','screened_activities')},indent=2))

if __name__=='__main__':main()
