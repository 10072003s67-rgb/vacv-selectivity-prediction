"""Fetch public CC50 assays in the same ChEMBL documents as VACV activities.

Matches by molecule and publication. Paired rows need cell-line/assay review before SI use.
"""
import csv
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1] / 'external'
API = 'https://www.ebi.ac.uk/chembl/api/data/activity.json'
CC_FIELDS = ('activity_id','molecule_chembl_id','canonical_smiles','document_chembl_id',
             'assay_chembl_id','assay_description','standard_relation','standard_value',
             'standard_units','data_validity_comment','potential_duplicate','source_url')
PAIR_FIELDS = ('molecule_chembl_id','document_chembl_id','vacv_activity_id','vacv_type',
               'vacv_value_nM','vacv_assay_chembl_id','vacv_assay_description',
               'cc50_activity_id','cc50_value_nM','cc50_assay_chembl_id',
               'cc50_assay_description','review_status')


def fetch(group):
    offset = 0
    while True:
        query = urllib.parse.urlencode({'document_chembl_id__in':','.join(group),
                                       'standard_type':'CC50','limit':1000,'offset':offset})
        for attempt in range(3):
            try:
                with urllib.request.urlopen(f'{API}?{query}',timeout=60) as response:
                    page = json.load(response)
                break
            except Exception:
                if attempt == 2: raise
                time.sleep(2)
        yield from page['activities']
        offset += len(page['activities'])
        if not page['page_meta']['next']: break


def write_csv(path, fields, records):
    with path.open('w',newline='',encoding='utf-8') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)


def main():
    vacv=list(csv.DictReader((ROOT/'chembl_vacv_activities.csv').open(encoding='utf-8')))
    documents=sorted({row['document_chembl_id'] for row in vacv if row['document_chembl_id']})
    cc={str(row['activity_id']):row for i in range(0,len(documents),20)
        for row in fetch(documents[i:i+20])}
    vacv_pairs={(row['document_chembl_id'],row['molecule_chembl_id']) for row in vacv}
    matched=sorted((row for row in cc.values() if
                    (row['document_chembl_id'],row['molecule_chembl_id']) in vacv_pairs),
                   key=lambda row:int(row['activity_id']))
    cc_rows=[]
    for row in matched:
        cc_rows.append({**{key:row.get(key) for key in CC_FIELDS if key!='source_url'},
            'source_url':f"https://www.ebi.ac.uk/chembl/compound_report_card/{row['molecule_chembl_id']}/"})
    write_csv(ROOT/'chembl_cc50_same_molecule_publication.csv',CC_FIELDS,cc_rows)
    exact={}
    for row in matched:
        if row['standard_relation']=='=' and row['standard_units']=='nM' and row['standard_value']:
            exact.setdefault((row['document_chembl_id'],row['molecule_chembl_id']),[]).append(row)
    candidates=[]
    for row in vacv:
        if (row['standard_type'] not in ('IC50','EC50') or row['standard_relation']!='='
                or row['standard_units']!='nM' or not row['standard_value']): continue
        for other in exact.get((row['document_chembl_id'],row['molecule_chembl_id']),[]):
            candidates.append({'molecule_chembl_id':row['molecule_chembl_id'],
                'document_chembl_id':row['document_chembl_id'],
                'vacv_activity_id':row['activity_id'],'vacv_type':row['standard_type'],
                'vacv_value_nM':row['standard_value'],'vacv_assay_chembl_id':row['assay_chembl_id'],
                'vacv_assay_description':row['assay_description'],
                'cc50_activity_id':other['activity_id'],'cc50_value_nM':other['standard_value'],
                'cc50_assay_chembl_id':other['assay_chembl_id'],
                'cc50_assay_description':other['assay_description'],
                'review_status':'review_cell_line_assay_and_leakage'})
    write_csv(ROOT/'chembl_cc50_vacv_pair_candidates.csv',PAIR_FIELDS,candidates)
    summary={'source':'ChEMBL activity API','vacv_documents':len(documents),
             'all_cc50_rows_in_vacv_publications':len(cc),
             'cc50_rows_same_molecule_and_publication':len(matched),
             'exact_nM_potency_cc50_pair_candidates':len(candidates),
             'unique_molecules_with_pair_candidates':len({x['molecule_chembl_id'] for x in candidates}),
             'accepted_SI_training_rows':0,
             'caution':'Pair candidates are not SI labels. Verify same cell line, assay conditions, units, structure and absence of train/test leakage before calculating SI.'}
    (ROOT/'chembl_cc50_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(summary)

if __name__=='__main__':main()
