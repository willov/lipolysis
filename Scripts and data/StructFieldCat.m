function S = StructFieldCat(S, T, dim)
if nargin<3 || isempty(dim)
    dim=1;
end
for i=1:size(T,2)
    fields = fieldnames(S);
    for k = 1:numel(fields)
        aField     = fields{k}; % EDIT: changed to {}
        if size(S.(aField),1) ==1 && size(S.(aField),2)>1 % If field is a horizontal 1D array
            S.(aField) = horzcat(S.(aField), T(i).(aField));
        elseif size(S.(aField),2) > 1 && size(S.(aField),1)>1 %if field is a vertical 1D array
            S.(aField) = vertcat(S.(aField), T(i).(aField));
        else
            S.(aField) = cat(dim, S.(aField), T(i).(aField));
        end
    end
end