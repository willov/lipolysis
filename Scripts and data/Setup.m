function []=Setup()

if ~exist('IQMsimulate','file')
    fprintf('\n\nFor these scripts to work, the IQM tools toolbox have to be compiled.\n')
    disp('To do this, a valid C-compiler is necessary')
    disp('If the compilation of the toolbox does not, it is likely that a valid C-compiler is missing.')
    choice = input('press enter to continue and compile the toolbox.');
    run('./IQM Tools/installIQMtoolsInitial.m')
end
%% Compile models
disp("Compiling models. Please wait.")
files = dir('*.txt');
for i = 1:length(files)
    IQMmakeMEXmodel(IQMmodel([files(i).folder '/' files(i).name]))
end
addpath(genpath('.'));
end
