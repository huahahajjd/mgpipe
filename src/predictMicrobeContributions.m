function [minFluxes,maxFluxes,fluxSpans,summary]=predictMicrobeContributions(modPath,varargin)
% Resumable per-sample microbial FVA. MATLAB R2019b compatible.
% Defaults: continuous rolling queue, 8 processes, IdleTimeout=Inf.
% Results are committed per sample, independently of batch completion.
% Use scheduling='rolling' to replace finished jobs immediately.
% metList contains VMH IDs; empty means all IEX reactions (original scope).
% returnTables=false avoids building large legacy cell tables in memory.
% See README.md for recovery import, memory admission and lean-model caching.
p=inputParser;positive=@(x) isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=1&&fix(x)==x;
p.addRequired('modPath',@ischar);
p.addParameter('resultsFolder',fullfile(pwd,'Contributions'),@ischar);
p.addParameter('cacheFolder','',@ischar);
p.addParameter('numWorkers',8,positive);p.addParameter('batchSize',8,positive);
p.addParameter('metList',{},@(x) iscellstr(x));
p.addParameter('targetRxns',{},@iscellstr);
p.addParameter('priorResultsFolder','',@ischar);
p.addParameter('numericalRetry',true,@islogical);
p.addParameter('scheduling','rolling',@(x) any(strcmp(x,{'batches','rolling'})));
p.addParameter('optPercentage',99.99,@(x) isnumeric(x)&&isscalar(x)&&x>0&&x<=100);
p.addParameter('reactionChunkSize',Inf,@(x) isequal(x,Inf)||positive(x));
p.addParameter('lpTimeLimit',Inf,@(x) isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('workerMemoryGiB',8,@(x) isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
p.addParameter('memoryBudgetGiB',0,@(x) isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('reserveMemoryGiB',8,@(x) isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0);
p.addParameter('returnTables',nargout>0,@islogical);p.addParameter('exportCSV',true,@islogical);
p.addParameter('recoveredFolder','',@ischar);p.addParameter('acceptRecovered',false,@islogical);
p.addParameter('dryRun',false,@islogical);p.addParameter('failOnError',true,@islogical);
p.addParameter('configureSolver',true,@islogical);p.addParameter('fvaFunction',@pmcFastFVA,@(x) isa(x,'function_handle'));
p.parse(modPath,varargin{:});o=p.Results;
assert(~isfinite(o.lpTimeLimit) || ~strcmp(func2str(o.fvaFunction),'fastFVA'), ...
    'PMC:UncheckedLimit','Finite LP limits require the status-checked pmcFastFVA backend.');
o.resultsFolder=absolute(o.resultsFolder);
if isempty(o.cacheFolder),o.cacheFolder=fullfile(o.resultsFolder,'lean_models');end
o.cacheFolder=absolute(o.cacheFolder);
for folder={o.resultsFolder,o.cacheFolder,fullfile(o.resultsFolder,'samples'),fullfile(o.resultsFolder,'logs')}
    if ~exist(folder{1},'dir'),mkdir(folder{1});end
end
% Prevent two coordinators writing the same results. Keep pool ownership separate.
lock=java.io.File(fullfile(o.resultsFolder,'.running.lock'));
assert(lock.mkdir(),'PMC:RunLocked','Results folder is locked by another run. Verify it has stopped before removing .running.lock.');
lockCleanup=onCleanup(@() lock.delete()); %#ok<NASGU>
files=dir(fullfile(modPath,'*.mat'));files=files(~[files.isdir]);
assert(~isempty(files),'PMC:NoModels','No MAT models in input folder.');
[~,order]=sort({files.name});files=files(order);
spec=struct('schema',2,'metList',{sort(unique(o.metList))},'targetRxns',{sort(unique(o.targetRxns))},'optPercentage',o.optPercentage, ...
    'backend',func2str(o.fvaFunction),'lpTimeLimit',o.lpTimeLimit);
analysisKey=pmcHash(jsonencode(spec));jobs=struct([]);done=false(numel(files),1);
for i=1:numel(files)
    source=absolute(fullfile(modPath,files(i).name));st=dir(source);
    sourceKey=pmcHash(sprintf('%s|%d|%.15g|lean-v1',source,st.bytes,st.datenum));
    [~,name]=fileparts(files(i).name);name=regexprep(name,'^\d+_','');
    job=struct('name',name,'source',source,'sourceKey',sourceKey, ...
        'key',pmcHash([sourceKey analysisKey]),'cache',fullfile(o.cacheFolder,[sourceKey '.mat']), ...
        'output',fullfile(o.resultsFolder,'samples',[name '.mat']));
    if i>1,assert(~any(strcmp({jobs(1:i-1).name},name)),'PMC:DuplicateModel','Duplicate normalized model name.');end
    if i==1,jobs=repmat(job,1,numel(files));else,jobs(i)=job;end
    if exist(job.output,'file')
        saved=load(job.output,'result');
        assert(isfield(saved,'result') && strcmp(saved.result.key,job.key),'PMC:CheckpointMismatch', ...
            'Existing result differs in source/settings: %s. Choose a new resultsFolder.',name);
        done(i)=strcmp(saved.result.state,'complete');
    else
        if ~isempty(o.priorResultsFolder),done(i)=pmcImportPrior(job,o);end
        if ~done(i) && o.acceptRecovered && ~isempty(o.recoveredFolder)
            done(i)=pmcImportRecovered(job,o);
        end
    end
end
budget=o.memoryBudgetGiB;
if budget==0
    budget=availableMemoryGiB()-o.reserveMemoryGiB;
end
active=min(o.numWorkers,floor(budget/o.workerMemoryGiB));
assert(active>=1,'PMC:Memory','Insufficient available memory for one worker with configured reserve.');
summary=struct('total',numel(jobs),'alreadyComplete',sum(done),'completed',sum(done), ...
    'numWorkers',o.numWorkers,'maxActive',active,'batchSize',o.batchSize,'scheduling',o.scheduling, ...
    'memoryBudgetGiB',budget,'workerMemoryGiB',o.workerMemoryGiB,'failed',{{}},'analysisKey',analysisKey);
fprintf('Models %d; resume %d; remaining %d; pool %d; active limit %d; scheduling %s.\n', ...
    numel(jobs),sum(done),sum(~done),o.numWorkers,active,o.scheduling);
if strcmp(o.scheduling,'rolling'),summary.batchSize=[];end
if active<o.numWorkers
    warning('PMC:MemoryAdmission','Memory budget limits active samples to %d; pool size remains %d. Adjust only after measuring peak worker RSS.',active,o.numWorkers);
end
minFluxes={};maxFluxes={};fluxSpans={};
pmcAtomicSave(fullfile(o.resultsFolder,'run_manifest.mat'),struct('jobs',jobs,'options',o,'summary',summary));
if o.dryRun,return;end
pending=find(~done);
if ~isempty(pending)
    pool=gcp('nocreate');
    if isempty(pool),pool=parpool('local',o.numWorkers,'IdleTimeout',Inf);
    else
        assert(pool.NumWorkers==o.numWorkers,'PMC:PoolSize','Existing pool has %d workers, requested %d. Resize explicitly before starting.',pool.NumWorkers,o.numWorkers);
        pool.IdleTimeout=Inf;
    end
    assert(isinf(pool.IdleTimeout),'PMC:PoolTimeout','Could not disable idle timeout.');
    environment=[];
    if o.configureSolver,environment=getEnvironment();end
    if strcmp(o.scheduling,'rolling'),groups={pending};
    else
        groups=arrayfun(@(first) pending(first:min(first+o.batchSize-1,end)), ...
            1:o.batchSize:numel(pending),'UniformOutput',false);
    end
    for b=1:numel(groups)
        assert(pool.Connected,'PMC:PoolLost','Pool disconnected; completed samples are on disk. No serial fallback.');
        if strcmp(o.scheduling,'rolling')
            fprintf('[%s] Start continuous queue: %d samples, up to %d active.\n',datestr(now,30),numel(groups{b}),active);
        else
            fprintf('[%s] Dispatch group %d/%d (%d samples).\n',datestr(now,30),b,numel(groups),numel(groups{b}));
        end
        summary=runGroup(pool,jobs(groups{b}),o,environment,active,summary);
        pmcAtomicSave(fullfile(o.resultsFolder,'run_summary.mat'),struct('summary',summary));
        if o.failOnError && ~isempty(summary.failed)
            error('PMC:SampleFailed','A sample failed. Completed samples remain saved. See logs; no unconstrained FBA fallback was used.');
        end
    end
end
if summary.completed==numel(jobs) && (o.exportCSV || o.returnTables)
    [minFluxes,maxFluxes,fluxSpans]=pmcExport(jobs,o);
end
pmcAtomicSave(fullfile(o.resultsFolder,'run_summary.mat'),struct('summary',summary));
end

function summary=runGroup(pool,jobs,o,environment,active,summary)
count=numel(jobs);futures=parallel.FevalFuture.empty;
holder=containers.Map({'futures'},{futures});
cancelCleanup=onCleanup(@() cancelHeld(holder)); %#ok<NASGU>
next=1;finished=0;slotJob=[];
while next<=min(active,count)
    futures(next)=parfeval(pool,@pmcRunSample,1,jobs(next),o,environment);holder('futures')=futures;slotJob(next)=next;next=next+1;
end
while finished<count
    [slot,state]=fetchNext(futures);
    jobIndex=slotJob(slot);finished=finished+1; %#ok<NASGU>
    if strcmp(state.state,'complete')
        summary.completed=summary.completed+1;
    else
        summary.failed{end+1}=state.name;
        fprintf(2,'%s\n',state.message);
    end
    pmcAtomicSave(fullfile(o.resultsFolder,'run_summary.mat'),struct('summary',summary));
    fprintf('[%s] %s %s (%.1fs), %d/%d dispatched samples finished.\n',datestr(now,30),state.name,state.state,state.seconds,finished,count);
    if next<=count
        futures(slot)=parfeval(pool,@pmcRunSample,1,jobs(next),o,environment);holder('futures')=futures;slotJob(slot)=next;next=next+1;
    end
end
end
function cancelHeld(holder)
futures=holder('futures');if ~isempty(futures),cancel(futures);end
end
function result=absolute(path)
result=char(java.io.File(path).getCanonicalPath());
end
function gib=availableMemoryGiB()
if isunix && exist('/proc/meminfo','file')
    token=regexp(fileread('/proc/meminfo'),'MemAvailable:\s+(\d+)','tokens','once');
    if ~isempty(token),gib=str2double(token{1})/2^20;return;end
end
error('PMC:MemoryBudget','Specify memoryBudgetGiB on systems without /proc/meminfo.');
end
