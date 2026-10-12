"""Check exact TIC identities in publisher-checksummed TARS default catalogue."""
import csv, datetime, hashlib, io, json, pathlib, zipfile
import requests
ROOT=pathlib.Path(__file__).resolve().parents[1]
OUT=ROOT/'data/tars-aggregate-review'
TICS={201308024,201308023,355767420,281593180,149392671,150104580}
def main():
    record=json.loads((OUT/'record.json').read_text())
    item=next(f for f in record['files'] if f['key']=='tars_table_2.csv.zip')
    dest=OUT/item['key']; md5=hashlib.md5(); sha=hashlib.sha256(); n=0
    if not dest.exists():
        with requests.get(item['links']['self'],stream=True,timeout=(30,180)) as response:
            response.raise_for_status()
            with dest.with_suffix('.part').open('wb') as handle:
                for block in response.iter_content(1024*1024):
                    handle.write(block)
        dest.with_suffix('.part').rename(dest)
    with dest.open('rb') as handle:
        for block in iter(lambda:handle.read(1024*1024),b''):
            md5.update(block); sha.update(block); n+=len(block)
    if 'md5:'+md5.hexdigest()!=item['checksum'] or n!=item['size']:
        raise ValueError('Publisher checksum/length mismatch')
    (OUT/'catalogue.provenance.json').write_text(json.dumps({'url':item['links']['self'],'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'bytes':n,'sha256':sha.hexdigest(),'publisher_checksum':item['checksum'],'publisher_checksum_verified':True},indent=2))
    matches=[]; gaia_matches=[]; count=0
    gaia_ids={"4906011953900670976","4906017863775668992"}
    with zipfile.ZipFile(dest) as archive:
        members=archive.namelist(); print(members,flush=True)
        names=[name for name in members if name.endswith('.csv')]
        if len(names)!=1: raise ValueError('Ambiguous CSV members')
        with archive.open(names[0]) as raw:
            reader=csv.DictReader(io.TextIOWrapper(raw,encoding='utf-8-sig'))
            columns=reader.fieldnames; print(columns,flush=True)
            ticcol=next(c for c in columns if c.lower() in ('tic','ticid','tic_id'))
            for row in reader:
                count+=1
                if int(row[ticcol]) in TICS: matches.append(row)
                if row["dr2_source_id"] in gaia_ids or row["dr3_source_id"] in gaia_ids: gaia_matches.append(row)
    result={'catalogue':'TARS table 2 default adopted rotation-period catalogue','row_count':count,'columns':columns,'exact_tics_queried':sorted(TICS),'matches':matches,'exact_gaia_ids_queried':sorted(gaia_ids),'exact_gaia_matches':gaia_matches,'source_sha256':hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),'unmatched_tics':sorted(TICS-{int(r[ticcol]) for r in matches}),'heldout_flux_opened':False,'scope':'Exact TIC membership only; nonmatches do not establish novelty or absence from full table 4 parent measurements.'}
    (OUT/'results.json').write_text(json.dumps(result,indent=2)); print(json.dumps(result,indent=2))
if __name__=='__main__': main()
