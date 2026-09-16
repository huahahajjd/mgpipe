function [lo,hi,obj,ret]=pmcChunkMock(~,~,~,~,rxns,varargin)
% Environment marker makes a repeatable failure between two checkpoints.
marker=fullfile(tempdir,'pmc_test_chunk_fail');
if any(contains(rxns,'_b[u]tr')) && exist(marker,'file'),error('PMC:InjectedChunk','Second chunk failure');end
if any(contains(rxns,'_a[u]tr')) && ~exist(marker,'file'),error('PMC:RepeatedChunk','Completed chunk was recalculated');end
lo=-ones(numel(rxns),1);hi=ones(numel(rxns),1);obj=1;ret=0;
end
