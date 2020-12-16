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
  run('./Scripts and data/Setup.m')
  choice = input('\n\n\nWhich figures do you want to plot? \n  1. Figure 2 and 3\n  2. Figure 4.\n  3. Figure 5, 6 and S1. \nChoice (1, 2 or 3): ');
  choice2 = input('\nWhat resolution do you want? \nHigh resolution (as in the papers) or low resolution? \nWarning: High resolution will take quite a long time to simulate. \n  1. High resolution\n  2. Low resolution\nChoice (1 or 2): ');
  
  if ~ismember(choice,[1 2 3]) && ~ismember(choice2,[1 2])
    error('Invalid choice of plots. Invalid choice of resolution')
  elseif ~ismember(choice,[1 2 3]) && ismember(choice2,[1 2])
    error('Invalid choice of plots. Valid choice of resolution')
  elseif ismember(choice,[1 2 3]) && ~ismember(choice2,[1 2])
    error('Valid choice of plots. Invalid choice of resolution')
  elseif choice2==1
    res=0.01;
  elseif choice2==2
    res=1;
  else
    error('Invalid choice.')
  end
  
  cd('Scripts and data')
  if choice == 1
    PlotUncertainty(0,0,res)
  elseif choice == 2
    PlotUncertainty(0,1, res)
  elseif choice == 3
    PlotUncertainty(1,0, res)
  else
    disp('Invalid choice if plots')
  end
  cd(currentDir)
catch err
  disp('Something went wrong.')
  cd(currentDir)
  rethrow(err)
end
