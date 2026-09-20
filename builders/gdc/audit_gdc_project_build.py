"""Independent audit of a built GDC project (run from the dir holding gdc_mirror/): re-fetch source files by a different code path and compare values, counts, column ids and clinical totals. Usage: python audit_gdc_project_build.py TCGA-CHOL [--mirror gdc_mirror]"""
import sys, json, io, gzip, random, subprocess, boto3, pandas as pd, numpy as np, pyarrow.parquet as pq
proj=sys.argv[1]; mirror=(sys.argv[sys.argv.index("--mirror")+1] if "--mirror" in sys.argv else "gdc_mirror"); out=f"{mirror}/{proj}"
rep=json.load(open(f"{out}/BUILD_REPORT.json"))
print("== report keys:", {k:v for k,v in rep.items() if k!='tables'})
tot=0
for name,t in rep['tables'].items(): tot+=t['bytes']; print(f"  {name:55s} {t['rows']:>9d} x {t['columns']:>5d}  {t['bytes']/1e6:8.2f} MB")
print("total MB", round(tot/1e6,1))
man=pd.read_parquet(f"{out}/MANIFEST.parquet")
s3=boto3.client('s3',region_name='us-east-1')
def fetch(r):
    b=s3.get_object(Bucket=r.s3_bucket,Key=r.s3_key)['Body'].read()
    return gzip.decompress(b) if r.file_name.endswith('.gz') else b
random.seed(3)
# --- RNA: 3 random files re-read directly, compare 5 random genes
tpm=pd.read_parquet(f"{out}/rna_seq/tpm_unstranded.parquet").set_index('gene_id')
cnt=pd.read_parquet(f"{out}/rna_seq/counts_unstranded.parquet").set_index('gene_id')
smp=pd.read_parquet(f"{out}/rna_seq/samples.parquet")
star=man[man.workflow_type=='STAR - Counts']
assert tpm.shape[1]==len(star)==cnt.shape[1], (tpm.shape,len(star))
assert tpm.shape[0]==60660, tpm.shape
bad=0
for _,r in star.sample(3,random_state=1).iterrows():
    src=pd.read_csv(io.BytesIO(fetch(r)),sep='\t',comment='#').iloc[4:].set_index('gene_id')
    col=smp.loc[smp.file_id==r.file_id,'column'].iloc[0]
    genes=random.sample(list(src.index),5)
    for g in genes:
        if not (np.isclose(src.loc[g,'tpm_unstranded'],tpm.loc[g,col]) and src.loc[g,'unstranded']==cnt.loc[g,col]): bad+=1
    print(f"  rna check {col}: {len(genes)} genes {'ok' if bad==0 else 'MISMATCH'}; dtype {tpm[col].dtype}/{cnt[col].dtype}")
# --- column ids are tumor aliquots? (TCGA barcode chars 14-15: 01-09 tumor, 10-19 normal)
for tbl in ['rna_seq/tpm_unstranded','copy_number/gene_level_ascat3','copy_number/gene_level_ascat2','mirna/rpm','protein/rppa_protein_expression']:
    try: cols=[c for c in pq.read_schema(f"{out}/{tbl}.parquet").names if c.startswith('TCGA-')]
    except Exception as e: print('  skip',tbl,e); continue
    codes=[c.split('-')[3][:2] for c in cols if len(c.split('-'))>3]
    normal=[c for c in codes if c.isdigit() and 10<=int(c)<=19]
    print(f"  {tbl}: {len(cols)} TCGA columns, sample codes {sorted(set(codes))}, normal-coded columns: {len(normal)}")
# --- MAF: row count equals sum of source rows for 2 files
maf=pq.read_metadata(f"{out}/somatic_mutations/masked_somatic_mutations.parquet")
mafman=man[man.data_type=='Masked Somatic Mutation']
mafpq=pd.read_parquet(f"{out}/somatic_mutations/masked_somatic_mutations.parquet",columns=['file_id','Hugo_Symbol','Variant_Classification','t_alt_count'])
for _,r in mafman.sample(2,random_state=2).iterrows():
    src=pd.read_csv(io.BytesIO(fetch(r)),sep='\t',comment='#',dtype=str)
    n=(mafpq.file_id==r.file_id).sum(); print(f"  maf {r.file_id[:8]}: source rows {len(src)} parquet rows {n} {'ok' if n==len(src) else 'MISMATCH'}")
print(f"  maf total rows {maf.num_rows}, cols {maf.num_columns}; t_alt_count dtype {mafpq.t_alt_count.dtype}; top classes {mafpq.Variant_Classification.value_counts().head(3).to_dict()}")
# --- clinical vs GDC case count
import urllib.request, urllib.parse
q=urllib.parse.urlencode({"filters":json.dumps({"op":"=","content":{"field":"project.project_id","value":proj}}),"size":0})
tot_cases=json.load(urllib.request.urlopen("https://api.gdc.cancer.gov/cases?"+q))['data']['pagination']['total']
cases=pd.read_parquet(f"{out}/clinical/cases.parquet"); dx=pd.read_parquet(f"{out}/clinical/diagnoses.parquet")
print(f"  clinical: cases {len(cases)} (GDC says {tot_cases}), diagnoses {len(dx)}, cols cases={len(cases.columns)} dx={len(dx.columns)}; vital_status {cases['demographic_vital_status'].value_counts().to_dict() if 'demographic_vital_status' in cases else 'MISSING'}")
print("  cases cols sample:", [c for c in cases.columns][:25])
print("  dx cols sample:", [c for c in dx.columns if 'stage' in c or 'diagnosis' in c][:10])
# --- samples table joins
sam=pd.read_parquet(f"{out}/clinical/aliquots.parquet")
missing=set(smp.aliquot_submitter_id)-set(sam.aliquot_submitter_id); print(f"  rna aliquots not in clinical/aliquots: {len(missing)}")
print("  dup/realign notes:", rep.get('duplicate_aliquots'), rep.get('realigned_files'))
