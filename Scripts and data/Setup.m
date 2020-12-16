function []=Setup()

if ~exist('IQMsimulate','file')
  choice = input('For these scripts to work, the IQM tools toolbox have to be compiled. \nThis needs a valid C-compiler.\nPress enter to continue.');
  run('./IQM Tools/installIQMtoolsInitial.m')
end
%% Compile models

addpath(genpath('.'));
end
