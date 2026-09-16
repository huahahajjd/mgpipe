function accepted=pmcImportPrior(job,o)
% Migrate a verified subset of complete results or completed reaction chunks.
accepted=false;
persistent folder old
if isempty(folder) || ~strcmp(folder,o.priorResultsFolder)
    old=load(fullfile(o.priorResultsFolder,'run_manifest.mat'),'jobs','options');
    assert(old.options.optPercentage==o.optPercentage,'PMC:PriorSettings','Different optimality fraction.');
    assert(isempty(old.options.metList),'PMC:PriorScope','Prior run must contain all IEX targets.');
    assert(~isfield(old.options,'targetRxns') || isempty(old.options.targetRxns),'PMC:PriorScope','Prior run is already restricted.');
    assert(strcmp(func2str(old.options.fvaFunction),'pmcFastFVA'),'PMC:PriorBackend','Prior results require checked backend.');
    folder=o.priorResultsFolder;
end
idx=find(strcmp({old.jobs.name},job.name));if isempty(idx),return;end
assert(isscalar(idx));prior=old.jobs(idx);
assert(strcmp(prior.sourceKey,job.sourceKey),'PMC:PriorSource','Source model changed since prior run.');
file=fullfile(folder,'samples',[job.name '.mat']);
if exist(file,'file')
    x=load(file,'result');r=x.result;
    assert(strcmp(r.key,prior.key) && strcmp(r.state,'complete'),'PMC:PriorResult','Invalid prior result.');
    selected=pmcSelectReactions(r.rxns,o.metList,o.targetRxns);
    [~,index]=ismember(selected,r.rxns);
    r.key=job.key;r.sourceKey=job.sourceKey;r.rxns=selected;
    r.minFlux=r.minFlux(index);r.maxFlux=r.maxFlux(index);r.migratedFrom=file;
    pmcAtomicSave(job.output,struct('result',r));accepted=true;
elseif exist([file '.progress.mat'],'file')
    p=load([file '.progress.mat']);
    assert(strcmp(p.key,prior.key),'PMC:PriorProgress','Invalid prior checkpoint key.');
    selected=pmcSelectReactions(p.rxns,o.metList,o.targetRxns);[~,index]=ismember(selected,p.rxns);
    minFlux=p.minFlux(index);maxFlux=p.maxFlux(index);completed=logical(p.completed(index));
    assert(all(isfinite(minFlux(completed))) && all(isfinite(maxFlux(completed))) && ...
        all(minFlux(completed)<=maxFlux(completed)+1e-6),'PMC:PriorValues','Invalid completed values.');
    if all(completed)
        r=struct('key',job.key,'sourceKey',job.sourceKey,'modelFile',job.name,'source',job.source, ...
            'rxns',{selected},'minFlux',minFlux,'maxFlux',maxFlux,'state','complete', ...
            'optPercentage',o.optPercentage,'origin','computed','solverReturnCode',0, ...
            'seconds',0,'schemaVersion',2,'migratedFrom',[file '.progress.mat']);
        pmcAtomicSave(job.output,struct('result',r));accepted=true;
    elseif ~exist([job.output '.progress.mat'],'file')
        pmcAtomicSave([job.output '.progress.mat'],struct('key',job.key,'rxns',{selected}, ...
            'minFlux',minFlux,'maxFlux',maxFlux,'completed',completed));
    end
end
end
