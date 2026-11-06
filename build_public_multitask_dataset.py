"""Build one public VACV multitask table; optional private competition files only filter overlaps.

Endpoint labels must be modelled separately. Derived SI labels remain provisional.
"""
import argparse,csv,json,re
from collections import defaultdict,Counter
from pathlib import Path
from rdkit import Chem,rdBase
from rdkit.Chem.MolStandardize import rdMolStandardize

ROOT=Path(__file__).resolve().parents[1]
EXT=ROOT/'external'
FIELDS=('record_id','molecule_chembl_id','canonical_smiles','inchikey','endpoint',
        'value','unit','label_quality','cell_line','document_chembl_id',
        'assay_chembl_id','cc50_assay_chembl_id','assay_description',
        'cc50_assay_description','source_url')
CELL_PATTERNS={'HFF':r'\bhff\b','BSC1':r'\bbsc[ -]?1\b','CV1':r'\bcv[ -]?1\b',
               'Vero':r'\bvero\b','MRC5':r'\bmrc[ -]?5\b','HeLa':r'\bhela\b',
               'BHK21':r'\bbhk[ -]?21\b','A549':r'\ba549\b'}

def cell(text):
    matches=[name for name,pattern in CELL_PATTERNS.items() if re.search(pattern,text or '',re.I)]
    return matches[0] if len(matches)==1 else ''

def structure(smiles):
    try:
        m=Chem.MolFromSmiles(smiles)
        if m is None:return None
        m=rdMolStandardize.FragmentParent(m)
        return Chem.MolToSmiles(m,canonical=True,isomericSmiles=True),Chem.MolToInchiKey(m)
    except Exception:return None

def excluded_keys(paths):
    keys=set()
    for path in paths:
        with Path(path).open(newline='',encoding='utf-8-sig') as stream:
            for row in csv.DictReader(stream):
                for column in ('canonical_smiles','SMILES','smiles_raw'):
                    if row.get(column):
                        parsed=structure(row[column])
                        if parsed:keys.add(parsed[1])
                        break
    return keys

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--exclude-smiles-csv',action='append',default=[],help='Private train/test CSV; used only for exclusion')
    args=parser.parse_args()
    excluded=excluded_keys(args.exclude_smiles_csv)
    vacv=list(csv.DictReader((EXT/'chembl_vacv_activities.csv').open(encoding='utf-8')))
    cc=list(csv.DictReader((EXT/'chembl_cc50_same_molecule_publication.csv').open(encoding='utf-8')))
    pairs=list(csv.DictReader((EXT/'chembl_cc50_vacv_pair_candidates.csv').open(encoding='utf-8')))
    structure_by_mol={}
    for row in vacv+cc:
        mid=row['molecule_chembl_id']
        if mid not in structure_by_mol and row.get('canonical_smiles'):
            structure_by_mol[mid]=structure(row['canonical_smiles'])
    output=[];seen=set()
    for row in vacv+cc:
        endpoint=row.get('standard_type','CC50')
        if endpoint not in ('IC50','EC50','CC50'):continue
        if row['standard_relation']!='=' or row['standard_units']!='nM':continue
        if row.get('data_validity_comment') or row.get('potential_duplicate') in ('1','True','true'):continue
        try:value=float(row['standard_value'])/1000
        except (TypeError,ValueError):continue
        if not 0<value<1e9:continue
        parsed=structure_by_mol.get(row['molecule_chembl_id'])
        if not parsed or parsed[1] in excluded:continue
        activity_id=row['activity_id']
        if activity_id in seen:continue
        seen.add(activity_id)
        output.append(dict(record_id='CHEMBL:'+activity_id,
            molecule_chembl_id=row['molecule_chembl_id'],canonical_smiles=parsed[0],
            inchikey=parsed[1],endpoint=endpoint,value=f'{value:.12g}',unit='uM',
            label_quality='measured_exact',cell_line=cell(row.get('assay_description')),
            document_chembl_id=row['document_chembl_id'],assay_chembl_id=row['assay_chembl_id'],
            cc50_assay_chembl_id='',assay_description=row.get('assay_description',''),
            cc50_assay_description='',source_url=row.get('source_url','')))
    matched=defaultdict(list)
    for row in pairs:
        if row['vacv_type']!='IC50':continue
        first=cell(row['vacv_assay_description']);second=cell(row['cc50_assay_description'])
        if not first or first!=second:continue
        matched[(row['document_chembl_id'],row['molecule_chembl_id'])].append((row,first))
    for rows in matched.values():
        if len(rows)!=1:continue
        row,cell_line=rows[0]
        parsed=structure_by_mol.get(row['molecule_chembl_id'])
        if not parsed or parsed[1] in excluded:continue
        try:si=float(row['cc50_value_nM'])/float(row['vacv_value_nM'])
        except (ValueError,ZeroDivisionError):continue
        if not 0<si<1e9:continue
        output.append(dict(record_id='DERIVED:'+row['vacv_activity_id']+':'+row['cc50_activity_id'],
            molecule_chembl_id=row['molecule_chembl_id'],canonical_smiles=parsed[0],
            inchikey=parsed[1],endpoint='SI',value=f'{si:.12g}',unit='ratio',
            label_quality='provisional_same_publication_cell_line',cell_line=cell_line,
            document_chembl_id=row['document_chembl_id'],assay_chembl_id=row['vacv_assay_chembl_id'],
            cc50_assay_chembl_id=row['cc50_assay_chembl_id'],
            assay_description=row['vacv_assay_description'],
            cc50_assay_description=row['cc50_assay_description'],
            source_url='https://www.ebi.ac.uk/chembl/compound_report_card/'+row['molecule_chembl_id']+'/'))
    output.sort(key=lambda x:(x['endpoint'],x['record_id']))
    path=EXT/'vacv_multitask_training.csv'
    with path.open('w',newline='',encoding='utf-8') as stream:
        writer=csv.DictWriter(stream,fieldnames=FIELDS);writer.writeheader();writer.writerows(output)
    summary={'rows':len(output),'unique_molecules':len({x['inchikey'] for x in output}),
             'endpoint_counts':dict(Counter(x['endpoint'] for x in output)),
             'excluded_competition_structures':len(excluded),
             'rdkit_version':rdBase.rdkitVersion,
             'warning':'SI rows are provisional derived labels. Model CC50, IC50, EC50 separately. Do not treat IC50/EC50/CC50 as SI or merge endpoint values. Competition train/test were used only to exclude molecule overlap and are not published.'}
    (EXT/'vacv_multitask_training_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(summary)

if __name__=='__main__':main()
