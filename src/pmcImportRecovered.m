function accepted=pmcImportRecovered(job,o)
% Recovery import is explicit because the original solver codes were not saved.
accepted=false;
assert(o.optPercentage==99.99 && any(strcmp(func2str(o.fvaFunction),{'fastFVA','pmcFastFVA'})), ...
    'PMC:RecoverySettings','Recovered outputs are for original fastFVA at 99.99%% only.');
persistent folder manifest
if isempty(folder) || ~strcmp(folder,o.recoveredFolder)
    report=jsondecode(fileread(fullfile(o.recoveredFolder,'matlab_validation.json')));
    assert(report.passed && report.numericSha256Verified,'Recovery validation not passed.');
    manifest=jsondecode(fileread(fullfile(o.recoveredFolder,'recovery_manifest.json')));
    folder=o.recoveredFolder;
end
names=regexprep({manifest.modelFile},'^\d+_','');
index=find(strcmp(names,[job.name '.mat']));
if isempty(index),return;end
assert(isscalar(index),'Duplicate recovered model.');record=manifest(index);
file=fullfile(o.recoveredFolder,'recovered_models',record.file);
fid=fopen(file,'r');assert(fid>=0,'Cannot read recovery file.');bytes=fread(fid,Inf,'*uint8');fclose(fid);
md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);
hash=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
assert(strcmp(hash,record.matSha256),'PMC:RecoveryHash','Recovered MAT checksum mismatch.');
r=load(file);
assert(strcmp(regexprep(r.modelFile,'^\d+_',''),[job.name '.mat']),'Recovered identity mismatch.');
rxns=pmcSelectReactions(r.rxns,o.metList,o.targetRxns);[found,idx]=ismember(rxns,r.rxns);assert(all(found));
result=struct('key',job.key,'sourceKey',job.sourceKey,'modelFile',job.name,'source',job.source, ...
    'rxns',{rxns},'minFlux',r.minFlux(idx),'maxFlux',r.maxFlux(idx),'state','complete', ...
    'optPercentage',99.99,'origin','recovered','solverReturnCode',NaN, ...
    'note','Byte-validated recovery; original solver codes unavailable; source identity accepted by explicit import.', ...
    'seconds',0,'schemaVersion',1);
pmcAtomicSave(job.output,struct('result',result));accepted=true;
end
