function test_pmc_installed()
here=fileparts(mfilename('fullpath'));root=fileparts(here);
pmcConfigureCobra('/mnt/nfs/wangchao/softwares/cobratoolbox', ...
    '/home/wangchao/workdir/softwares/cplex/cplex129/CPLEX_Studio');
addpath(here,'-begin');
target='/mnt/nfs/wangchao/softwares/cobratoolbox/src/analysis/multiSpecies/microbiomeModelingToolbox/additionalAnalysis';
addpath(target,'-begin');
assert(strcmp(which('predictMicrobeContributions'),fullfile(target,'predictMicrobeContributions.m')));
assert(strcmp(which('pmcFastFVA'),fullfile(target,'pmcFastFVA.m')));
fprintf('Installed function: %s\n',which('predictMicrobeContributions'));
test_pmc_export();
folder=tempname;mkdir(folder);cleanup=onCleanup(@() rmdir(folder,'s')); %#ok<NASGU>
model=struct('S',sparse([1,-1]),'rxns',{{'in';'bug_IEX_ac[u]tr'}},'mets',{{'m'}}, ...
    'lb',[0;0],'ub',[10;10],'c',[0;1],'b',0,'csense','E', ...
    'C',sparse([0,1]),'d',4,'dsense','L','osenseStr','max');
control=struct('PARALLELMODE',1,'THREADS',1,'AUXROOTTHREADS',0);
[lo,hi,~,ret]=pmcFastFVA(model,99.99,'max','ibm_cplex',model.rxns(2),'A',control);
assert(ret==0 && abs(lo-3.9996)<1e-7 && abs(hi-4)<1e-7);
% Compare original 12.8 and new 12.9 kernels on one worker (no nested pool).
environment=getEnvironment();f=parfeval(gcp(),@pmcCompareOriginal,1,model,environment);
agreement=fetchOutputs(f);assert(agreement<1e-7);
models=fullfile(folder,'models');mkdir(models);
for i=1:8,save(fullfile(models,sprintf('sample_%d.mat',i)),'model');end
[~,~,~,summary]=predictMicrobeContributions(models,'resultsFolder',fullfile(folder,'results'), ...
    'numWorkers',8,'batchSize',8,'workerMemoryGiB',1,'memoryBudgetGiB',128, ...
    'returnTables',false,'exportCSV',false,'reactionChunkSize',1);
assert(summary.completed==8);
% Reject a deliberately iteration-limited LP instead of reporting it complete.
n=100;limited=model;limited.S=sparse(1:n,1:n,1,n,n+1)+sparse(1:n,2:n+1,-1,n,n+1);
limited=rmfield(limited,{'C','d','dsense'});limited.lb=zeros(n+1,1);limited.ub=ones(n+1,1);
limited.c=zeros(n+1,1);limited.b=zeros(n,1);limited.csense=repmat('E',n,1);
limited.rxns=arrayfun(@(i) sprintf('r%d',i),(1:n+1)','UniformOutput',false);
control.ITLIM=1;control.PREIND=0;control.LPMETHOD=1;
[~,~,~,limitedStatus]=pmcFastFVA(limited,99.99,'max','ibm_cplex',limited.rxns(end),'A',control);
assert(limitedStatus~=0,'Iteration-limited solution incorrectly accepted.');
fprintf('PASS: analytical coupled LP, original-kernel agreement %.3g, 8 real CPLEX workers, limited solve rejected (code %d).\n',agreement,limitedStatus);
end
