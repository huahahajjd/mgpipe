function [model, info] = pmcReadModel(filename)
% Read LP fields only. For v7.3, follow HDF5 references without loading annotations.
t = tic;
fid=fopen(filename,'r'); assert(fid>=0,'Cannot open model: %s',filename);
header=char(fread(fid,128,'*uint8')'); fclose(fid);
fields=pmcModelFields(); model=struct();
if contains(header,'MATLAB 7.3 MAT-file')
    file=H5F.open(filename,'H5F_ACC_RDONLY','H5P_DEFAULT');
    cleanup=onCleanup(@() H5F.close(file));
    root=H5G.open(file,'/'); rootCleanup=onCleanup(@() H5G.close(root));
    names=groupNames(root);
    containers={};
    for i=1:numel(names)
        if startsWith(names{i},'#'), continue; end
        obj=H5O.open(root,names{i},'H5P_DEFAULT');
        isModel=isGroup(obj) && hasChild(obj,'S') && hasChild(obj,'rxns');
        H5O.close(obj);
        if isModel, containers{end+1}=names{i}; end %#ok<AGROW>
    end
    if all(ismember({'S','rxns'},names)), group=root;
    else
        assert(numel(containers)==1,'PMC:ModelVariable','Expected exactly one model in %s.',filename);
        group=H5G.open(root,containers{1});
        groupCleanup=onCleanup(@() H5G.close(group)); %#ok<NASGU>
        names=groupNames(group);
    end
    for i=1:numel(fields)
        name=fields{i};
        if ismember(name,names)
            obj=H5O.open(group,name,'H5P_DEFAULT');
            objCleanup=onCleanup(@() H5O.close(obj));
            model.(name)=readObject(obj);
            clear objCleanup
        end
    end
    info.loader='selected HDF5 fields';
else
    variables=whos('-file',filename);
    if all(ismember({'S','rxns'},{variables.name}))
        model=load(filename,fields{:});
    else
        candidates={variables(strcmp({variables.class},'struct')).name};
        models={};
        for i=1:numel(candidates)
            s=load(filename,candidates{i}); x=s.(candidates{i});
            if isscalar(x) && isfield(x,'S') && isfield(x,'rxns'), models{end+1}=x; end %#ok<AGROW>
        end
        assert(numel(models)==1,'PMC:ModelVariable','Expected exactly one model in %s.',filename);
        model=models{1};
        model=rmfield(model,setdiff(fieldnames(model),fields));
    end
    info.loader='MAT load then retain LP fields';
end
assert(isfield(model,'S') && isfield(model,'rxns'),'PMC:MissingModel','Missing S/rxns.');
model.rxns=model.rxns(:);
assert(size(model.S,2)==numel(model.rxns),'PMC:Dimensions','S/rxns mismatch.');
assert(numel(unique(model.rxns))==numel(model.rxns),'PMC:DuplicateRxns','Duplicate reaction IDs.');
if isfield(model,'F') && nnz(model.F)>0
    error('PMC:QuadraticModel','Nonzero quadratic objective is not supported by LP FVA.');
end
if isfield(model,'vartype') && any(model.vartype~='C')
    error('PMC:IntegerModel','Integer variables are not supported by LP FVA.');
end
% The installed fastFVA indexes results by reaction count, not extra variables.
if isfield(model,'E') && size(model.E,2)>0
    error('PMC:ExtraVariables','This fastFVA adapter cannot safely handle extra variables E.');
end
if ~isfield(model,'b'),model.b=zeros(size(model.S,1),1);end
if ~isfield(model,'c'),model.c=zeros(size(model.S,2),1);end
if ~isfield(model,'csense'),model.csense=repmat('E',size(model.S,1),1);end
for name={'lb','ub','b','c','csense','d','dsense','dxdt'}
    if isfield(model,name{1}),v=model.(name{1});model.(name{1})=v(:);end
end
assert(isfield(model,'lb') && isfield(model,'ub'),'Missing flux bounds.');
assert(numel(model.lb)==size(model.S,2) && numel(model.ub)==size(model.S,2) && numel(model.c)==size(model.S,2), ...
    'PMC:Dimensions','Bounds/objective size mismatch.');
assert(numel(model.b)==size(model.S,1) && numel(model.csense)==size(model.S,1), ...
    'PMC:Dimensions','Mass-balance right hand side/sense size mismatch.');
if isfield(model,'C')
    assert(isfield(model,'d') && isfield(model,'dsense'),'PMC:Coupling','C requires d and dsense; never discard coupling constraints.');
    assert(size(model.C,2)==size(model.S,2) && size(model.C,1)==numel(model.d) && numel(model.d)==numel(model.dsense));
end
w=whos('model'); info.modelBytes=w.bytes;info.loadSeconds=toc(t);
info.reactions=numel(model.rxns);info.metabolites=size(model.S,1);
info.couplingRows=0;if isfield(model,'C'),info.couplingRows=size(model.C,1);end
end

function names=groupNames(group)
n=H5G.get_num_objs(group);names=cell(1,n);
for i=1:n,names{i}=H5G.get_objname_by_idx(group,i-1);end
end
function yes=hasChild(group,name)
yes=H5L.exists(group,name,'H5P_DEFAULT');
end
function v=attribute(obj,name)
a=H5A.open_name(obj,name);cleanup=onCleanup(@() H5A.close(a)); %#ok<NASGU>
v=H5A.read(a);
end
function value=readObject(obj)
attributes=attributeNames(obj);
if isGroup(obj)
    assert(ismember('MATLAB_sparse',attributes),'PMC:HDFType','Unsupported group in LP fields.');
    nr=double(attribute(obj,'MATLAB_sparse'));
    ir=readChild(obj,'ir');jc=readChild(obj,'jc');
    if hasChild(obj,'data'),v=readChild(obj,'data');else,v=[];end
    nc=numel(jc)-1;
    cols=repelem((1:nc)',diff(double(jc(:))));
    value=sparse(double(ir(:))+1,cols,double(v(:)),nr,nc);
    return
end
cls='';if ismember('MATLAB_class',attributes),cls=char(attribute(obj,'MATLAB_class')');cls=cls(:)';end
raw=H5D.read(obj);
if ismember('MATLAB_empty',attributes) && attribute(obj,'MATLAB_empty')
    dims=double(raw(:)');
    if strcmp(cls,'cell'),value=cell(dims);elseif strcmp(cls,'char'),value=char(zeros(dims));else,value=zeros(dims,cls);end
    return
end
if strcmp(cls,'cell')
    space=H5D.get_space(obj);[~,dims]=H5S.get_simple_extent_dims(space);H5S.close(space);
    dims=fliplr(double(dims));if numel(dims)==1,dims=[dims,1];end
    refs=reshape(raw,8,[]);value=cell(dims);
    for k=1:numel(value)
        target=H5R.dereference(obj,'H5R_OBJECT',refs(:,k));
        cleanup=onCleanup(@() H5O.close(target));
        value{k}=readObject(target);clear cleanup
    end
elseif strcmp(cls,'char'),value=char(raw);
elseif strcmp(cls,'logical'),value=logical(raw);
else
    assert(isnumeric(raw),'PMC:HDFType','Unsupported numeric representation.');value=raw;
end
end
function v=readChild(group,name)
d=H5D.open(group,name);cleanup=onCleanup(@() H5D.close(d));v=H5D.read(d); %#ok<NASGU>
end
function yes=isGroup(obj)
kind=H5I.get_type(obj);
yes=isequal(kind,H5ML.get_constant_value('H5I_GROUP')) || strcmp(kind,'H5I_GROUP');
end
function names=attributeNames(obj)
n=H5A.get_num_attrs(obj);names=cell(1,n);
for i=1:n
    a=H5A.open_idx(obj,i-1);names{i}=H5A.get_name(a);H5A.close(a);
end
end
