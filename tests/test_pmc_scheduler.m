function test_pmc_scheduler()
here=fileparts(mfilename('fullpath'));addpath(here,fileparts(here));
root=tempname;mkdir(root);cleanup=onCleanup(@() rmdir(root,'s')); %#ok<NASGU>
models=fullfile(root,'models');mkdir(models);
model=struct('S',sparse([1,-1]),'rxns',{{'in';'bug_IEX_ac[u]tr'}},'mets',{{'m'}}, ...
    'lb',[0;0],'ub',[10;10],'c',[0;1],'b',0,'csense','E');
for i=1:10,save(fullfile(models,sprintf('model_%02d.mat',i)),'model');end
opts={'scheduling','batches','numWorkers',8,'batchSize',8,'configureSolver',false,'fvaFunction',@pmcMockFVA, ...
    'memoryBudgetGiB',128,'workerMemoryGiB',1,'returnTables',false,'exportCSV',true};
out=fullfile(root,'results');
[~,~,~,s]=predictMicrobeContributions(models,'resultsFolder',out,opts{:});
assert(s.completed==10);pool=gcp('nocreate');assert(pool.NumWorkers==8 && isinf(pool.IdleTimeout));
originalPool=pool;
[~,~,~,s]=predictMicrobeContributions(models,'resultsFolder',out,opts{:});
assert(s.alreadyComplete==10 && isequal(originalPool,gcp('nocreate')));
numeric=load(fullfile(out,'fluxes_numeric.mat'));assert(all(full(numeric.minFluxes)==-1));
% Model edits cannot silently reuse an old result.
model.ub(1)=11;save(fullfile(models,'model_01.mat'),'model');
failed=false;try,predictMicrobeContributions(models,'resultsFolder',out,opts{:});catch e,failed=strcmp(e.identifier,'PMC:CheckpointMismatch');end
assert(failed);
% Continuous scheduling, partial final group, per-sample failure and retry.
model.ub(1)=12345;save(fullfile(models,'model_02.mat'),'model');
out2=fullfile(root,'rolling');
[~,~,~,s]=predictMicrobeContributions(models,'resultsFolder',out2,opts{:},'scheduling','rolling','failOnError',false);
assert(s.completed==9 && numel(s.failed)==1);
assert(~exist(fullfile(out2,'samples','model_02.mat'),'file'));
assert(isequal(originalPool,gcp('nocreate')));
fprintf('PASS: 8-worker persistent pool, batches 8+2, rolling dispatch, resume, invalidation and failure isolation.\n');
end
