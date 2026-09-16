function fields = pmcModelFields()
% Keep every field used to define the LP, including community coupling rows.
fields = {'S','C','d','dsense','b','c','lb','ub','csense','rxns','mets', ...
    'E','D','evarlb','evarub','evarc','dxdt','osense','osenseStr','F','vartype'};
end
