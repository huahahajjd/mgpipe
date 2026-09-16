"""Back up v2 and atomically install the validated v3 MATLAB source files."""
from pathlib import Path
import datetime,hashlib,json,os,shutil
base=Path('/mnt/nfs/wangchao/my_course/untreated_MS/replace_data/step4_modeling/predict_microbe_contributions')
stage=base/'refactor_v3_20260915'
target=Path('/mnt/nfs/wangchao/softwares/cobratoolbox/src/analysis/multiSpecies/microbiomeModelingToolbox/additionalAnalysis')
old=json.loads((base/'refactor_20260914'/'installation_manifest.json').read_text())
for row in old['files']:
 p=Path(row['path']);assert hashlib.sha256(p.read_bytes()).hexdigest()==row['sha256'],'Unexpected installed edit: '+str(p)
assert (base/'reference_targets_20260915.mat').is_file()
backup=base/'code_backups'/datetime.datetime.now().strftime('before_v3_%Y%m%d_%H%M%S');backup.mkdir(parents=True)
names=['predictMicrobeContributions.m','pmcSelectReactions.m','pmcRunSample.m','pmcImportRecovered.m','pmcImportPrior.m','prepare_reference_targets.m']
pairs=[(stage/n,target/n) for n in names]+[(stage/'run_predict_microbe_contributions_v3.m',base/'run_predict_microbe_contributions_v3.m'),(stage/'README_v3.md',base/'README_predictMicrobeContributions_v3.md')]
records=[]
for src,dst in pairs:
 assert src.is_file() and src.stat().st_size
 record={'path':str(dst),'sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'backup':None}
 if dst.exists():
  saved=backup/dst.name;shutil.copy2(str(dst),str(saved));record['backup']=str(saved)
  record['previousSha256']=hashlib.sha256(saved.read_bytes()).hexdigest()
 tmp=dst.with_name(dst.name+'.installing');shutil.copy2(str(src),str(tmp));os.replace(str(tmp),str(dst))
 assert hashlib.sha256(dst.read_bytes()).hexdigest()==record['sha256'];records.append(record)
manifest={'installedAt':datetime.datetime.now().isoformat(),'backupDirectory':str(backup),'files':records}
(stage/'installation_v3.json').write_text(json.dumps(manifest,indent=2));(backup/'installation_v3.json').write_text(json.dumps(manifest,indent=2))
print(json.dumps({'installed':len(records),'backup':str(backup)},indent=2))
