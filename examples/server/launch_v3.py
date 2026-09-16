"""Start the authorized rolling production job detached from SSH."""
from pathlib import Path
import datetime,hashlib,json,subprocess,sys
base=Path('/mnt/nfs/wangchao/my_course/untreated_MS/replace_data/step4_modeling/predict_microbe_contributions')
results=base/'run_reference_rolling_20260915'
assert not (results/'.running.lock').exists(),'Existing run lock; inspect before starting another task.'
manifest=json.loads((base/'refactor_v3_20260915'/'installation_v3.json').read_text())
for r in manifest['files']:assert hashlib.sha256(Path(r['path']).read_bytes()).hexdigest()==r['sha256'],r['path']
launch=base/'production_launches'/datetime.datetime.now().strftime('reference_rolling_%Y%m%d_%H%M%S');launch.mkdir(parents=True)
expression="addpath('%s'); summary=run_predict_microbe_contributions_v3(); fprintf('PMC_QUEUE_FINISHED completed=%%d total=%%d failed=%%d\\n',summary.completed,summary.total,numel(summary.failed)); if ~isempty(summary.failed), error('PMC:Incomplete','Queue finished; unresolved samples remain. See run_summary.mat.'); end; disp('PMC_PRODUCTION_COMPLETE');"%base
command=['/mnt/nfs/wangchao/softwares/R2019b_matlab/bin/matlab','-nodisplay','-nosplash','-batch',expression]
runner=launch/'runner.py'
runner.write_text('''import datetime,fcntl,json,subprocess
from pathlib import Path
base=Path(%r);launch=Path(%r);results=Path(%r)
lock=open(str(results/'.launcher.lock'),'a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
command=%r
p=subprocess.Popen(command,cwd=str(base),stdin=subprocess.DEVNULL)
record={'matlabPid':p.pid,'startedAt':datetime.datetime.now().isoformat(),'command':command,'log':str(launch/'matlab.log')}
(launch/'process.json').write_text(json.dumps(record,indent=2))
code=p.wait();record.update(exitCode=code,finishedAt=datetime.datetime.now().isoformat())
(launch/'exit_status.json').write_text(json.dumps(record,indent=2))
'''%(str(base),str(launch),str(results),command))
with open(str(launch/'matlab.log'),'ab',buffering=0) as log:
 p=subprocess.Popen([sys.executable,str(runner)],cwd=str(base),stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
record={'runnerPid':p.pid,'launchDirectory':str(launch),'resultsDirectory':str(results),'version':'v3-reference-rolling','reference':'ERR3472679','referenceTargetCount':7928}
(launch/'launch.json').write_text(json.dumps(record,indent=2));(base/'latest_production_launch.json').write_text(json.dumps(record,indent=2));print(json.dumps(record,indent=2))
