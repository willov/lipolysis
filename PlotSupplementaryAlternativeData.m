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
    run('./Scripts and data/Setup.m')
    res = 0.01;
    cd('Scripts and data')
    PlotPL('lipolysis', 0)
    PlotUncertainty(0,res, 'lipolysis', 0)
    cd(currentDir)
catch err
    disp('Something went wrong.')
    cd(currentDir)
    rethrow(err)
end