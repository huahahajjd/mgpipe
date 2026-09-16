function test_pmc_export()
root=tempname;mkdir(root);cleanup=onCleanup(@() rmdir(root,'s')); %#ok<NASGU>
for j=1:2
    result=struct('rxns',{{'panBug_IEX_ac[u]tr'}},'minFlux',(-1)^j,'maxFlux',2); %#ok<NASGU>
    jobs(j)=struct('output',fullfile(root,sprintf('%d.mat',j)),'name',sprintf('microbiota_model_samp_s%d',j)); %#ok<AGROW>
    save(jobs(j).output,'result');
end
o=struct('resultsFolder',root,'returnTables',true,'exportCSV',true);
[lo,hi,span]=pmcExport(jobs,o);
assert(size(lo,1)==2 && strcmp(lo{2,2},'-1') && strcmp(lo{2,3},'1'));
assert(ischar(hi{2,2}) && isnumeric(span{2,2}) && span{2,2}==3);
n=load(fullfile(root,'fluxes_numeric.mat'));assert(isequal(full(n.minFluxes),[-1,1]));
assert(exist(fullfile(root,'Microbe_Flux_Spans.csv'),'file')==2);
fprintf('PASS: opposite signs cannot delete a nonzero row; legacy output types and numeric MAT/CSV export.\n');
end
