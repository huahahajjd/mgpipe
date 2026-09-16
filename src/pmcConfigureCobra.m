function pmcConfigureCobra(cobraRoot,cplexRoot)
% Configure the existing installation without a network/update initialization.
addpath(genpath(cobraRoot));
addpath(fullfile(cplexRoot,'cplex','matlab','x86-64_linux'),'-begin');
setenv('ILOG_CPLEX_PATH',cplexRoot);setenv('CPLEX_STUDIO_DIR129',cplexRoot);setenv('CPLEX_STUDIO_DIR',cplexRoot);
global SOLVERS OPT_PROB_TYPES CBTDIR ENV_VARS ILOG_CPLEX_PATH
global CBT_LP_SOLVER CBT_MILP_SOLVER CBT_QP_SOLVER CBT_MIQP_SOLVER CBT_NLP_SOLVER CBT_EP_SOLVER CBT_CLP_SOLVER
CBTDIR=cobraRoot;ILOG_CPLEX_PATH=cplexRoot;
OPT_PROB_TYPES={'LP','MILP','QP','MIQP','NLP','EP','CLP'};
SOLVERS=struct('ibm_cplex',struct('type',{{'LP','MILP','QP','MIQP'}},'categ','legacy','installed',true,'working',true));
ENV_VARS=struct('STATUS',1,'printLevel',false);
CBT_LP_SOLVER='ibm_cplex';CBT_MILP_SOLVER='';CBT_QP_SOLVER='';CBT_MIQP_SOLVER='';
CBT_NLP_SOLVER='';CBT_EP_SOLVER='';CBT_CLP_SOLVER='';
assert(changeCobraSolver('ibm_cplex','LP',1,0),'CPLEX initialization failed.');
end
