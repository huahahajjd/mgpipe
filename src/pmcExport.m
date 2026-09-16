function [minFluxes,maxFluxes,fluxSpans]=pmcExport(jobs,o)
% Build numeric sparse tables once, after calculation. No growing string cells.
allIDs=cell(numel(jobs),1);names=cell(1,numel(jobs));
for j=1:numel(jobs)
    x=load(jobs(j).output,'result');allIDs{j}=x.result.rxns;
    names{j}=regexprep(jobs(j).name,'^microbiota_model_(samp|diet)_','');
end
[ids,~,row]=unique(vertcat(allIDs{:}),'stable');
counts=cellfun(@numel,allIDs);total=sum(counts);cols=zeros(total,1);lo=cols;hi=cols;offset=0;
for j=1:numel(jobs)
    x=load(jobs(j).output,'result');index=offset+(1:counts(j));
    cols(index)=j;lo(index)=x.result.minFlux;hi(index)=x.result.maxFlux;offset=offset+counts(j);
end
low=sparse(row,cols,lo,numel(ids),numel(jobs));high=sparse(row,cols,hi,numel(ids),numel(jobs));
spans=high-low;
labels=regexprep(strrep(strrep(ids,'_IEX',''),'[u]tr',''),'^pan','');
keepLow=any(abs(low)>1e-7,2);keepHigh=any(abs(high)>1e-7,2);keepSpan=any(abs(spans)>1e-7,2);
if o.exportCSV
    writeCSV(fullfile(o.resultsFolder,'Microbe_Secretion.csv'),labels(keepLow),names,low(keepLow,:));
    writeCSV(fullfile(o.resultsFolder,'Microbe_Uptake.csv'),labels(keepHigh),names,high(keepHigh,:));
    writeCSV(fullfile(o.resultsFolder,'Microbe_Flux_Spans.csv'),labels(keepSpan),names,spans(keepSpan,:));
end
pmcAtomicSave(fullfile(o.resultsFolder,'fluxes_numeric.mat'),struct('rxns',{ids},'samples',{names}, ...
    'minFluxes',low,'maxFluxes',high,'fluxSpans',spans));
minFluxes={};maxFluxes={};fluxSpans={};
if o.returnTables
    cells=(sum(keepLow)+sum(keepHigh)+sum(keepSpan))*numel(names);
    assert(cells<=1e7,'PMC:LegacyTablesTooLarge', ...
        'Legacy cell tables would exceed 10 million cells. Numeric MAT and CSV are saved; use returnTables=false.');
    minFluxes=legacy(labels(keepLow),names,low(keepLow,:),true);
    maxFluxes=legacy(labels(keepHigh),names,high(keepHigh,:),true);
    fluxSpans=legacy(labels(keepSpan),names,spans(keepSpan,:),false);
end
end
function value=legacy(ids,names,matrix,asText)
values=num2cell(full(matrix));
if asText,values=cellfun(@(v) sprintf('%.17g',v),values,'UniformOutput',false);end
value=[{''},names;ids,values];
end
function writeCSV(file,ids,names,matrix)
tmp=[file '.partial'];fid=fopen(tmp,'w');assert(fid>=0,'Cannot write CSV.');
cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'reaction');for j=1:numel(names),fprintf(fid,',"%s"',strrep(names{j},'"','""'));end
fprintf(fid,'\n');
for first=1:512:numel(ids)
    last=min(first+511,numel(ids));block=full(matrix(first:last,:));
    for r=1:size(block,1)
        fprintf(fid,'"%s"',strrep(ids{first+r-1},'"','""'));
        fprintf(fid,',%.17g',block(r,:));fprintf(fid,'\n');
    end
end
clear cleanup
[ok,msg]=movefile(tmp,file,'f');assert(ok,'%s',msg);
end
