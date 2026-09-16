function rxns=pmcSelectReactions(allRxns,metList,targetRxns)
% Exact VMH suffix matching avoids selecting similarly named metabolites.
if isempty(metList)
    mask=contains(allRxns,'IEX_');
else
    mask=false(size(allRxns));
    for i=1:numel(metList)
        suffix=['IEX_' metList{i} '[u]tr'];
        mask=mask | strcmp(allRxns,suffix) | endsWith(allRxns,['_' suffix]);
    end
end
if nargin>=3 && ~isempty(targetRxns),mask=mask & ismember(allRxns,targetRxns);end
rxns=allRxns(mask);rxns=rxns(:);
end
