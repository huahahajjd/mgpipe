function test_pmc_reference()
% Exercise exact strain-specific targets and migration of completed subsets.
root=tempname;mkdir(root);cleanup=onCleanup(@() rmdir(root,'s')); %#ok<NASGU>
rxns={'in';'bug_IEX_ac[u]tr';'bug_IEX_but[u]tr'};
selected=pmcSelectReactions([rxns;{'other_IEX_ac[u]tr'}],{},rxns(2));assert(isequal(selected,rxns(2)));
model=struct('S',sparse([1,-1,-1]),'rxns',{rxns},'mets',{{'m'}},'lb',[0;0;0], ...
    'ub',[10;10;10],'c',[0;1;1],'b',0,'csense','E','C',sparse([0,1,1]),'d',4,'dsense','L');
models=fullfile(root,'models');mkdir(models);for i=1:2,save(fullfile(models,sprintf('s%d.mat',i)),'model');end
old=fullfile(root,'old');cache=fullfile(root,'cache');
opts={'cacheFolder',cache,'memoryBudgetGiB',128,'returnTables',false,'exportCSV',false,'reactionChunkSize',1};
[~,~,~,s]=predictMicrobeContributions(models,'resultsFolder',old,opts{:});assert(s.completed==2);
f=fullfile(old,'samples','s2.mat');z=load(f,'result');r=z.result;
pmcAtomicSave([f '.progress.mat'],struct('key',r.key,'rxns',{r.rxns},'minFlux',r.minFlux, ...
    'maxFlux',r.maxFlux,'completed',[true;false]));delete(f);
out=fullfile(root,'selected');
[~,~,~,s]=predictMicrobeContributions(models,'resultsFolder',out,'priorResultsFolder',old,'targetRxns',rxns(2),opts{:});
assert(s.alreadyComplete==2 && s.completed==2 && isempty(s.batchSize));
z=load(fullfile(out,'samples','s2.mat'));assert(isequal(z.result.rxns,rxns(2)) && z.result.maxFlux==r.maxFlux(1));
% Select the uncomputed second reaction: migrate its flag and solve only it.
[~,~,~,s]=predictMicrobeContributions(models,'resultsFolder',fullfile(root,'pending'), ...
    'priorResultsFolder',old,'targetRxns',rxns(3),opts{:});assert(s.alreadyComplete==1 && s.completed==2);
% Verify numerical-stability parameters are recognized by the actual backend.
control=struct('THREADS',1,'NUMERICALEMPHASIS',1,'SCAIND',1,'ADVIND',1);
[lo,hi,~,ret]=pmcFastFVA(model,99.99,'max','ibm_cplex',rxns(2),'A',control);
assert(ret==0 && abs(lo)<1e-6 && abs(hi-4)<1e-6);
fprintf('PASS: exact strain-specific targets, complete/partial migration, pending subset solve, numerical controls.\n');
end
