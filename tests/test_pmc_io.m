function test_pmc_io()
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
folder=tempname;mkdir(folder);cleanup=onCleanup(@() rmdir(folder,'s')); %#ok<NASGU>
model=struct('S',sparse([1,-1,0;0,1,-1]),'rxns',{{'a';'bug_IEX_ac[u]tr';'bug_IEX_acald[u]tr'}}, ...
    'mets',{{'m1';'m2'}},'lb',[0;0;0],'ub',[10;10;10],'c',[0;0;1], ...
    'b',[0;0],'csense','EE','C',sparse([0,1,-2]),'d',0,'dsense','L','osenseStr','max', ...
    'rxnNames',{{'not needed';'not needed';'not needed'}});
for version={'-v7','-v7.3'}
    file=fullfile(folder,['model' version{1} '.mat']);save(file,'model',version{1});
    [lean,info]=pmcReadModel(file);
    for f=setdiff(fieldnames(model),{'rxnNames'})'
        a=model.(f{1});b=lean.(f{1});
        if ischar(a) && isvector(a) && ~strcmp(f{1},'osenseStr'),a=a(:);end
        assert(isequaln(a,b),['Mismatch in ' f{1}]);
    end
    assert(~isfield(lean,'rxnNames'));disp(info);
end
assert(isequal(pmcSelectReactions(model.rxns,{'ac'}),{'bug_IEX_ac[u]tr'}));
assert(isempty(pmcSelectReactions(model.rxns,{'missing'})));
file=fullfile(folder,'atomic.mat');pmcAtomicSave(file,struct('value',rand(3)));
s=load(file);assert(isequal(size(s.value),[3,3]));
fprintf('PASS: selected-field MAT readers, coupling preservation, exact targets, atomic save.\n');
end
