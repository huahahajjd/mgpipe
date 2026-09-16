function [lo,hi,obj,ret]=pmcMockFVA(model,~,~,~,rxns,varargin)
% Deterministic backend for scheduling/checkpoint tests; never used in production.
if model.ub(1)==12345,error('PMC:InjectedFailure','Intentional failure');end
pause(0.1);lo=-ones(numel(rxns),1);hi=2*ones(numel(rxns),1);obj=1;ret=0;
end
