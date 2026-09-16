function [lo,hi,objective,ret]=pmcFastFVA(model,percentage,sense,solver,rxns,matrixAS,control,varargin)
% Small-output adapter to the status-checked fastFVA kernel.
% No nested parallel pool, no repeated version discovery or directory changes.
assert(exist('pmcCplexFVA','file')==3,'PMC:BackendMissing','Run build_pmc_backend once on the server.');
assert(strcmp(solver,'ibm_cplex') && strcmp(sense,'max'),'Only CPLEX/max is supported.');
assert(~isempty(rxns),'PMC:EmptyTargets','Empty targets must be skipped before calling the backend.');
[found,index]=ismember(rxns,model.rxns);assert(all(found) && numel(unique(index))==numel(index));
problem=buildOptProblemFromModel(model);
% Match the original predictMicrobeContributions/fastFVA matrix selection.
if strcmp(matrixAS,'A') && isfield(model,'C')
    A=problem.A;b=problem.b;csense=problem.csense;c=problem.c;lb=problem.lb;ub=problem.ub;
else
    A=model.S;b=model.b;csense=model.csense;c=model.c;lb=model.lb;ub=model.ub;
end
assert(size(A,2)==numel(model.rxns),'PMC:ExtraVariables','Unsupported additional variables.');
if ~issparse(A),A=sparse(A);end
fields=fieldnames(control);values=zeros(numel(fields),1);
for i=1:numel(fields),values(i)=control.(fields{i});end
[allLo,allHi,objective,ret]=pmcCplexFVA(full(c(:)),A,full(b(:)),csense(:),full(lb(:)),full(ub(:)),percentage,-1, ...
    double(index(:)),1,control,values,2*ones(1,numel(index)),tempdir,0);
lo=allLo(index);hi=allHi(index);
end
