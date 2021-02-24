function [s] = StructsFieldCat(structs, dim)
%STRUCTSFIELDCAT Summary of this function goes here
%   Detailed explanation goes here
if nargin<2 || isempty(dim)
    dim=1;
end
s = structs(1);
for i = 2:length(structs)
    s=StructFieldCat(s, structs(i), dim);
end
end

