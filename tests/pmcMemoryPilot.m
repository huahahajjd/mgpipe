function report=pmcMemoryPilot(cacheFile,environment)
% Bounded resource pilot; a time-limited solve is diagnostic, never a result.
restoreEnvironment(environment);changeCobraSolver('ibm_cplex','LP',0,-1);
report=struct('pid',feature('getpid'),'before',fileread('/proc/self/status'));
t=tic;s=load(cacheFile,'model');report.loadSeconds=toc(t);model=s.model;clear s
rxns=pmcSelectReactions(model.rxns,{});report.target=rxns{1};
control=struct('PARALLELMODE',1,'THREADS',1,'AUXROOTTHREADS',0,'TILIM',20);
t=tic;[lo,hi,obj,ret]=pmcFastFVA(model,99.99,'max','ibm_cplex',rxns(1),'A',control);
report.solveSeconds=toc(t);report.ret=ret;report.min=lo;report.max=hi;report.objective=obj;
report.after=fileread('/proc/self/status');
end
