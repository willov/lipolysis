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
    run('./Scripts and data/Setup.m')
    fprintf('\n\nWhich figures do you want to plot?\n')
    disp('1. Parameter identifiability (Fig. 2)')
    disp('2. All figures for model agreement before diabetes was included (Fig. 3, 4, S1)')
    disp('3. Effects of removing insulin actions (Fig. 5)')
    disp('4. All figure after diabetes was introduced (Fig. 6, 7, S3, S4)')
    disp('5. Best agreement when insulin action 3 was removed (Fig. S2)')
    disp('Note that figures numbered 5X corresponds to supplementary figures SX.')
    choice = input('\nChoice (1, 2, 3, 4 or 5): ');
    
    if ~ismember(choice, [1,3,5])
        fprintf('\n\nWhat resolution do you want?\n')
        disp('High resolution (as in the papers) or low resolution (faster)?')
        disp('1. High resolution')
        disp('2. Low resolution')
        choice2 = input('\nChoice (1 or 2): ');
    else
        choice2 = 1;
    end
    fprintf('\n\n')
    if ~ismember(choice,[1 2 3 4 5]) && ~ismember(choice2,[1 2])
        error('Invalid choice of plots. Invalid choice of resolution')
    elseif ~ismember(choice,[1 2 3 4 5]) && ismember(choice2,[1 2])
        error('Invalid choice of plots. Valid choice of resolution')
    elseif ismember(choice,[1 2 3 4 5]) && ~ismember(choice2,[1 2])
        error('Valid choice of plots. Invalid choice of resolution')
    elseif choice2==1
        res=0.01;
    elseif choice2==2
        res=1;
    else
        error('Invalid choice.')
    end
    
    cd('Scripts and data')
    if choice ==1
        PlotPL
    elseif choice == 2
        PlotUncertainty(0,res)
    elseif choice == 3
        PlotInsulinActions
    elseif choice == 4
        PlotUncertainty(1, res)
    elseif choice == 5
        PlotUncertainty(0,res, 'lipolysis_noIns3')
    else
        disp('Invalid choice if plots')
    end
    cd(currentDir)
catch err
    disp('Something went wrong.')
    cd(currentDir)
    rethrow(err)
end
