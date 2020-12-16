function S = StructFieldCat(S, T, dim)
if nargin<3 || isempty(dim)
    dim=1;
end
for i=1:size(T,2)
    fields = fieldnames(S);
    for k = 1:numel(fields)
        aField     = fields{k}; 
        S.(aField) = cat(dim, S.(aField), T(i).(aField));
    end
end