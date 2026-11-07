"""Expand existing table with independent, exact ChEMBL measurements discovered via lek."""
import argparse, hashlib, json, math, re, time
from functools import lru_cache
from pathlib import Path
from urllib.request import urlopen
from urllib.parse import urlencode
import pandas as pd
from rdkit import Chem, rdBase
from rdkit.Chem.MolStandardize import rdMolStandardize

LEK='https://raw.githubusercontent.com/lenyagrigorenko/lek/f2c063fd70c5c427d09c6ef415bdeb38ed1320e6/data/external/vacv_external.csv'
BASE='https://raw.githubusercontent.com/10072003s67-rgb/vacv-selectivity-prediction/ae1eb43b48a147a538bc82b8359381b959735b24/vacv_training.csv'
API='https://www.ebi.ac.uk/chembl/api/data/activity.json'

@lru_cache(maxsize=None)
def structure(s):
    m=Chem.MolFromSmiles(str(s))
    if m is None: raise ValueError('Invalid SMILES')
    m=rdMolStandardize.FragmentParent(m)
    normalized=rdMolStandardize.Uncharger().uncharge(m)
    return Chem.MolToSmiles(m),Chem.MolToInchiKey(normalized)[:14]

def download(url,path):
    if not path.exists():
        for attempt in range(4):
            try:
                with urlopen(url,timeout=120) as response: body=response.read()
                path.write_bytes(body);break
            except Exception:
                if attempt==3: raise
                time.sleep(2**attempt)
    return path

def usable(a):
    return (a.get('standard_relation')=='=' and a.get('standard_type') in ('IC50','EC50','CC50')
        and a.get('standard_units') in ('nM','uM','M') and not a.get('data_validity_comment')
        and str(a.get('potential_duplicate')).lower() not in ('1','true'))

def main():
    p=argparse.ArgumentParser();p.add_argument('--train',required=True);p.add_argument('--test',required=True)
    p.add_argument('--cache',default='.cache/expansion');p.add_argument('--output',default='.')
    args=p.parse_args();cache=Path(args.cache);cache.mkdir(parents=True,exist_ok=True)
    out=Path(args.output);(out/'reports').mkdir(parents=True,exist_ok=True)
    basepath=download(BASE,cache/'base.csv');srcpath=download(LEK,cache/'lek.csv')
    base=pd.read_csv(basepath);src=pd.read_csv(srcpath,dtype=str).fillna('')
    excluded={structure(s)[1] for path in (args.train,args.test) for s in pd.read_csv(path).SMILES}
    seen=set(base.source_id); candidates={}
    for _,r in src.iterrows():
        if structure(r.smiles)[1] in excluded: continue
        for prefix in ('ic50','cc50'):
            if r[prefix+'_type'] not in ('IC50','EC50','CC50'):continue
            # Fetch each original activity: aggregated values in lek are never used as measurements.
            for aid in r[prefix+'_activity_ids'].split(','):
                if aid.isdigit() and 'CHEMBL:'+aid not in seen:
                    candidates.setdefault(aid,r.to_dict())
    activities=[];ids=sorted(candidates,key=int)
    for start in range(0,len(ids),50):
        batch=ids[start:start+50];url=API+'?'+urlencode({'activity_id__in':','.join(batch),'limit':1000})
        path=download(url,cache/('activities_'+hashlib.sha256(','.join(batch).encode()).hexdigest()[:16]+'.json'))
        data=json.loads(path.read_text());activities.extend(data['activities'])
    additions=[];provenance=[]
    for a in activities:
        sid='CHEMBL:'+str(a['activity_id']);r=candidates[str(a['activity_id'])]
        if sid in seen or not usable(a):continue
        if a.get('molecule_chembl_id')!=r['record_id'] or a.get('document_chembl_id')!=r['document_id']:continue
        if not a.get('canonical_smiles'):continue
        smiles,key=structure(a['canonical_smiles'])
        if key in excluded:continue
        endpoint=a['standard_type']
        if endpoint in ('IC50','EC50') and a.get('assay_chembl_id')!=r['assay_id']:continue
        if endpoint=='CC50' and a.get('assay_chembl_id')!=r['cc50_assay_id']:continue
        try:value=float(a['standard_value'])*{'nM':.001,'uM':1,'M':1e6}[a['standard_units']]
        except (ValueError,TypeError):continue
        if not math.isfinite(value) or value<=0:continue
        row=dict.fromkeys(base.columns,'')
        column='CC50_uM' if endpoint=='CC50' else endpoint+'_uM_VACV'
        row.update(Title=a['molecule_chembl_id'],SMILES=smiles,source_id=sid,label_quality='measured_exact')
        row[column]=value
        culture=r['cc50_cell_line'] if endpoint=='CC50' else r['cell_line']
        row['Culture']='' if culture in ('other','unknown') else culture
        additions.append(row);seen.add(sid)
        provenance.append(dict(source_id=sid,doi=r['doi'],document_id=a['document_chembl_id'],assay_id=a['assay_chembl_id'],endpoint=endpoint,source_url='https://www.ebi.ac.uk/chembl/api/data/activity/'+str(a['activity_id'])+'.json',discovery_source=LEK))
    added=pd.DataFrame(additions,columns=base.columns)
    result=pd.concat([base,added],ignore_index=True)
    assert len(result)==len(base)+len(added)
    assert not set(added.SMILES.map(lambda s:structure(s)[1]))&excluded
    assert added['SI VACV'].fillna('').eq('').all()
    # The published training table contains only SI labels and permitted inputs.
    # CC50 and IC50 are used solely while auditing external measurements.
    training_columns=['Title','strain','Culture','SMILES','SI VACV','source_id','label_quality']
    training=result.loc[result['SI VACV'].notna() & result['SI VACV'].astype(str).str.strip().ne(''),training_columns].copy()
    training.to_csv(out/'vacv_training.csv',index=False)
    pd.DataFrame(provenance).to_csv(out/'reports/added_measurements.csv',index=False)
    oldkeys={structure(s)[1] for s in base.SMILES}
    newkeys={structure(s)[1] for s in added.SMILES}
    report={'before_rows':len(base),'added_measurements':len(added),'after_rows':len(result),'training_rows':len(training),'training_columns':training_columns,'added_by_endpoint':pd.Series([r['endpoint'] for r in provenance]).value_counts().to_dict(),'new_molecules_connectivity':len(newkeys-oldkeys),'added_molecules_connectivity':len(newkeys),'new_si_labels':0,'excluded_competition_connectivities':len(excluded),'rdkit':rdBase.rdkitVersion,'base_url':BASE,'discovery_url':LEK,'base_sha256':hashlib.sha256(basepath.read_bytes()).hexdigest(),'discovery_sha256':hashlib.sha256(srcpath.read_bytes()).hexdigest(),'note':'The 50 measurements are retained for provenance only; they add no SI labels to the published training table. CC50 and IC50 are not training inputs.'}
    (out/'reports/expansion.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
