function pmcAtomicSave(filename,payload)
% One writer per sample; rename only after save and structural verification.
folder=fileparts(filename);if ~exist(folder,'dir'),mkdir(folder);end
tmp=[tempname(folder) '.mat'];cleanup=onCleanup(@() removeTemp(tmp));
% Small per-model payloads are much faster/smaller in v7 (especially cell IDs).
% Upgrade individual variables before reaching the 2 GiB v7 format limit.
fields=fieldnames(payload);version='-v7';
for i=1:numel(fields)
    value=payload.(fields{i});w=whos('value');
    if w.bytes>=1.8*2^30,version='-v7.3';break;end
end
save(tmp,'-struct','payload',version);
v=whos('-file',tmp);assert(numel(v)==numel(fieldnames(payload)),'PMC:Save','Incomplete MAT file.');
[ok,msg]=movefile(tmp,filename,'f');assert(ok,'PMC:Save','%s',msg);
end
function removeTemp(path)
if exist(path,'file'),delete(path);end
end
