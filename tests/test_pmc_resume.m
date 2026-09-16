function test_pmc_resume()
here=fileparts(mfilename('fullpath'));addpath(here,fileparts(here));
root=tempname;mkdir(root);cleanup=onCleanup(@() rmdir(root,'s')); %#ok<NASGU>
models=fullfile(root,'models');mkdir(models);
model=struct('S',sparse([1,-1,0]),'rxns',{{'in';'bug_IEX_a[u]tr';'bug_IEX_b[u]tr'}}, ...
    'mets',{{'m'}},'lb',[0;0;0],'ub',[10;10;10],'c',[0;1;0],'b',0,'csense','E');
save(fullfile(models,'model.mat'),'model');out=fullfile(root,'results');
marker=fullfile(tempdir,'pmc_test_chunk_fail');fid=fopen(marker,'w');fclose(fid);
markerCleanup=onCleanup(@() removeMarker(marker)); %#ok<NASGU>
opts={'configureSolver',false,'fvaFunction',@pmcChunkMock,'memoryBudgetGiB',128, ...
    'reactionChunkSize',1,'returnTables',false,'exportCSV',false,'failOnError',false};
[~,~,~,s]=predictMicrobeContributions(models,'resultsFolder',out,opts{:});assert(s.completed==0);
f=fullfile(out,'samples','model.mat.progress.mat');x=load(f);assert(isequal(x.completed,[true;false]));
delete(marker);
[~,~,~,s]=predictMicrobeContributions(models,'resultsFolder',out,opts{:});assert(s.completed==1 && ~exist(f,'file'));
% An absent target must create a valid empty result without invoking the solver.
[~,~,~,s]=predictMicrobeContributions(models,'resultsFolder',fullfile(root,'empty'),opts{:},'metList',{'absent'});
assert(s.completed==1);x=load(fullfile(root,'empty','samples','model.mat'));assert(isempty(x.result.rxns));
fprintf('PASS: resume skips completed reaction chunks, missing target never invokes solver.\n');
end
function removeMarker(file)
if exist(file,'file'),delete(file);end
end
