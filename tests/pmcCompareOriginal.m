function difference=pmcCompareOriginal(model,environment)
restoreEnvironment(environment);changeCobraSolver('ibm_cplex','LP',0,-1);
control=struct('PARALLELMODE',1,'THREADS',1,'AUXROOTTHREADS',0);
old=pwd;cleanup=onCleanup(@() cd(old)); %#ok<NASGU>
[lo1,hi1,~,ret1]=fastFVA(model,99.99,'max','ibm_cplex',model.rxns(2),'A',control,0,[],0);
[lo2,hi2,~,ret2]=pmcFastFVA(model,99.99,'max','ibm_cplex',model.rxns(2),'A',control);
assert(ret1==0 && ret2==0);difference=max(abs([lo1(:)-lo2(:);hi1(:)-hi2(:)]));
end
