% This file allow you to generate all the figures from the paper.
% For it to work, please change the current directory in MATLAB to the directory where the
% files are saved.
% If this file is run from the wrong directory, it will attempt to switch
% to the correct one, with permission from the user.


if ~(exist('./Model equations.txt','file') && exist('./Scripts and data','dir'))
    cdChoice=input('Warning! You are likely running the file from the wrong directory.\nDo you want to attempt to switch to the correct directory automatically?\n  1. yes\n  2. no\nChoice (1 or 2): ','s');
    if cdChoice == 1
        cd(fileparts(mfilename('fullpath')))
    elseif cdChoice == 2
        disp('Warning! The script will most likely not run from this directory.')
        disp('If you encounter errors, please change the directory to the folder with the files.' )
    else
        disp('Invalid choice')
    end
end
assert(exist('./Model equations.txt','file') && exist('./Scripts and data','dir'),'Error! Likely running from the wrong directory.')
currentDir = pwd;

try
    close all
    fprintf('\nThis version plots all the figures at the original resolution.')
    disp('This will take a while.')
    disp('The figures will be finished in batches, but all figures will be available when the script is done.')
    disp('Note that figures numbered 5X corresponds to supplementary figures SX.')
    choice = input('\npress enter to continue.');
    
    run('./Scripts and data/Setup.m')
    
    res = 0.01;
    
    cd('Scripts and data')
    fprintf('\n\nRecreating the parameter identifiability figure  (Fig. 2)\n')
    PlotPL
    fprintf('\n\nRecreating all plots related to figures before the addition of diabetes (Fig. 3, 4, S1)\n')
    PlotUncertainty(0,res)
    fprintf('\n\nRecreating the effects of removing insulin actions (Fig. 5)\n')
    PlotInsulinActions
    fprintf('\n\nRecreating all figure after diabetes was introduced (Fig. 6, 7, S3, S4)\n')
    PlotUncertainty(1, res)
    fprintf('\n\nRecreating the best agreement when insulin action 3 was removed (Fig. S2)\n')
    PlotUncertainty(0,res, 'lipolysis_noIns3')
    
    cd(currentDir)
catch err
    disp('Something went wrong.')
    cd(currentDir)
    rethrow(err)
end
