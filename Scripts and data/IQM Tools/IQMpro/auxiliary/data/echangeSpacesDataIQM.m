function [dataNew] = echangeSpacesDataIQM(data)
% This function removes all spaces in string elements in a matlab table
% and exchanges them agains ":::". This is sometimes needed when
% outdated software (like NONMEM) is unable to handle spaces in strings. I
% know this is a feature and not bug in the view of NONMEM users - I
% personally think this is not acceptable.
%
% The function restoreSpacesDataIQM() will reverse this change.
%
% [SYNTAX]
% [dataNew] = echangeSpacesDataIQM(data)
%
% [INPUT]
% data:             Dataset in MATLAB table format or string 
%
% [OUTPUT]
% dataNew:          Dataset in MATLAB table format or string with all spaces (' ')
%                   exchanged with ':::' 

% <<<COPYRIGHTSTATEMENT - IQM TOOLS PRO>>>

% Check input arguments
if istable(data),
    % Exchange all spaces in string variables with ':::' in dataset
    varNames = data.Properties.VariableNames;
    for k=1:length(varNames),
        if ~isnumeric(data.(varNames{k})),
            data.(varNames{k}) = strrep(data.(varNames{k}),' ',':::');
        end
    end
    % Assign output variable
    dataNew = data;
elseif ischar(data),
    dataNew = strrep(data,' ',':::');    
else
    error('Input argument is not a MATLAB table or a string.');
end



