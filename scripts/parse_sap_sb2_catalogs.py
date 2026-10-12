"""Reviewed reparse of preserved Traven SB2 queries; no network access."""
import io
from pathlib import Path
import warnings
from astropy.io.votable import parse
from astropy.io.votable.exceptions import W03, W50
from adastra.archive import write_json,utc_now,sha256_file
from review_qlp_catalogs import primitive

def main():
    out=Path('data/binary-prior-201308024')
    dest=out/'reparsed-results.json'
    if dest.exists():raise ValueError('Preserve previous parse')
    paths=[out/(name+'.xml') for name in ('galah-spectrum','galah-gaia')]
    write_json(out/'parse-amendment.json',{'created_utc':utc_now(),
        'source_sha256':sha256_file(__file__), 'inputs':{str(p):sha256_file(p) for p in paths},
        'reviewed_metadata_warning':"W50 log(cm.s**-2) corresponds to logarithmic surface gravity [cm/s2] in authoritative ReadMe; ID-generation warnings affect XML IDs, not field names. Values are not unit-converted.",
        'readme_sha256':sha256_file(out/'galah-sb2-ReadMe.txt'),
        'original_failures_preserved':True})
    results=[]
    for path in paths:
        with warnings.catch_warnings(record=True) as records:
            warnings.simplefilter('always');d=parse(io.BytesIO(path.read_bytes()),verify='warn')
        for w in records:
            msg=str(w.message)
            if not ((w.category is W03 and 'Implicitly generating an ID' in msg) or
                    (w.category is W50 and "Invalid unit string 'log(cm.s**-2)'" in msg)):
                raise ValueError('Unexpected warning: '+msg)
        infos=list(d.infos)+[i for resource in d.resources for i in resource.infos]
        if any((i.name or '').lower()=='error' or (i.value or '').upper() in ('ERROR','OVERFLOW') for i in infos):
            raise ValueError('Service error/overflow')
        tables=[]
        for table in d.iter_tables():
            if len(table.array)>=100:raise ValueError('Reached cap')
            tables.append({'catalogue':table.name,'rows':[
                {n:primitive(row[n]) for n in table.array.dtype.names} for row in table.array]})
        if 'J/A+A/638/A145/catalog' not in [t['catalogue'] for t in tables]:
            raise ValueError('Requested catalogue absent')
        results.append({'path':str(path),'status':'complete','tables':tables,
                        'warnings':[str(w.message) for w in records]})
        print(path.name,sum(len(t['rows']) for t in tables))
    write_json(dest,{'completed_utc':utc_now(),'queries':results})

if __name__=='__main__':main()
