function reference=prepare_reference_targets(censusFile,referenceFile,cacheFolder)
% Verify the smallest candidate against the actual model, then freeze its IDs.
census=jsondecode(fileread(censusFile));
assert(numel(census)==418,'PMC:CensusIncomplete','Expected all 418 source models.');
if isstruct(census),census=num2cell(census);end
n=cellfun(@(r) r.reactions,census);c=cellfun(@(r) r.couplingRows,census);
idx=find(n==min(n) & c==min(c),1);
assert(~isempty(idx),'PMC:NoJointMinimum','No sample simultaneously minimizes reactions and coupling rows; choose a criterion explicitly.');
record=census{idx};source=char(java.io.File(record.file).getCanonicalPath());
st=dir(source);assert(st.bytes==record.bytes,'PMC:SourceChanged','Census source size changed.');
[model,info]=pmcReadModel(source);
assert(numel(model.rxns)==record.reactions && info.couplingRows==record.couplingRows,'PMC:CensusMismatch','Metadata reader disagrees with model.');
targetRxns=pmcSelectReactions(model.rxns,{});assert(~isempty(targetRxns),'PMC:EmptyReference','Reference has no IEX reactions.');
reference=struct('source',source,'name',record.name,'reactions',record.reactions, ...
    'couplingRows',record.couplingRows,'targetRxns',{targetRxns},'targetCount',numel(targetRxns), ...
    'selection','Simultaneous minimum of reaction count and coupling-row count', ...
    'matching','Exact full reaction ID, including strain prefix', ...
    'targetHash',pmcHash(jsonencode(sort(targetRxns))));
pmcAtomicSave(referenceFile,struct('reference',reference));
sourceKey=pmcHash(sprintf('%s|%d|%.15g|lean-v1',source,st.bytes,st.datenum));
if ~exist(cacheFolder,'dir'),mkdir(cacheFolder);end
pmcAtomicSave(fullfile(cacheFolder,[sourceKey '.mat']),struct('model',model,'info',info,'sourceKey',sourceKey));
disp(rmfield(reference,'targetRxns'));
end
