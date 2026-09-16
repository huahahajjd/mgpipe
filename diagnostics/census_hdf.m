function records=census_hdf()
base='/*/predict_microbe_contributions';
files=dir(fullfile(base,'run_20260507_204155','ordered_models','*.mat'));records={};
for i=1:numel(files)
    path=fullfile(files(i).folder,files(i).name);fid=fopen(path,'r');head=char(fread(fid,128,'*uint8')');fclose(fid);
    if ~contains(head,'MATLAB 7.3'),continue;end
    f=H5F.open(path,'H5F_ACC_RDONLY','H5P_DEFAULT');cleanup=onCleanup(@() H5F.close(f));
    root=H5G.open(f,'/');rootCleanup=onCleanup(@() H5G.close(root));
    group=root;ownGroup=false;
    if ~H5L.exists(root,'S','H5P_DEFAULT')
        for j=0:H5G.get_num_objs(root)-1
            name=H5G.get_objname_by_idx(root,j);if startsWith(name,'#'),continue;end
            obj=H5O.open(root,name,'H5P_DEFAULT');kind=H5I.get_type(obj);
            isGroup=isequal(kind,H5ML.get_constant_value('H5I_GROUP')) || strcmp(kind,'H5I_GROUP');
            found=isGroup && H5L.exists(obj,'S','H5P_DEFAULT');H5O.close(obj);
            if found,group=H5G.open(root,name);ownGroup=true;break;end
        end
        assert(ownGroup,'No model');
    end
    [rows,cols]=shape(group,'S');coupling=0;
    if H5L.exists(group,'C','H5P_DEFAULT'),[coupling,~]=shape(group,'C');end
    if ownGroup,H5G.close(group);end
    records{end+1}=struct('name',files(i).name,'file',path,'bytes',files(i).bytes, ...
        'reactions',cols,'metaboliteRows',rows,'couplingRows',coupling); %#ok<AGROW>
    fprintf('HDF %d: %s, %d reactions, %d coupling rows\n',numel(records),files(i).name,cols,coupling);
    clear rootCleanup cleanup
end
assert(numel(records)==66);
fid=fopen(fullfile(base,'refactor_20260914','model_dimensions_v73.json'),'w');fprintf(fid,'%s',jsonencode(records));fclose(fid);
end
function [rows,cols]=shape(group,name)
g=H5G.open(group,name);cleanup=onCleanup(@() H5G.close(g)); %#ok<NASGU>
a=H5A.open_name(g,'MATLAB_sparse');rows=double(H5A.read(a));H5A.close(a);
d=H5D.open(g,'jc');space=H5D.get_space(d);[~,dims]=H5S.get_simple_extent_dims(space);cols=double(prod(dims))-1;
H5S.close(space);H5D.close(d);
end
