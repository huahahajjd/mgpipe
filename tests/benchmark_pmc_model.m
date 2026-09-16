function benchmark_pmc_model()
here=fileparts(mfilename('fullpath'));addpath(fileparts(here));
file='/home/wangchao/workdir/my_course/untreated_MS/replace_data/step4_modeling/part2_modeling/personlized_modeling/microbiota_model_samp_ERR6996065.mat';
[model,info]=pmcReadModel(file);w=whos('model');
fields=fieldnames(model);sizes=struct();
for i=1:numel(fields),v=model.(fields{i});z=whos('v');sizes.(fields{i})=z.bytes;end
out=fullfile(fileparts(here),'benchmark_largest_lean.mat');
t=tic;pmcAtomicSave(out,struct('model',model,'info',info));saveSeconds=toc(t);
t=tic;s=load(out,'model');reloadSeconds=toc(t);
assert(isequaln(model,s.model));
source=dir(file);target=dir(out);
report=struct('sourceBytes',source.bytes,'leanMemoryBytes',w.bytes,'cacheBytes',target.bytes, ...
    'loadSeconds',info.loadSeconds,'saveSeconds',saveSeconds,'reloadSeconds',reloadSeconds, ...
    'fieldBytes',sizes,'reactions',info.reactions,'couplingRows',info.couplingRows);
fid=fopen(fullfile(fileparts(here),'benchmark_report.json'),'w');fprintf(fid,'%s',jsonencode(report));fclose(fid);
disp(report);
end
