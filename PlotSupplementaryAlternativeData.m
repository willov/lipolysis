fprintf('\n\nWhat resolution do you want?\n')
disp('High resolution (as in the papers) or low resolution (faster)?')
disp('1. High resolution')
disp('2. Low resolution')
choice2 = input('\nChoice (1 or 2): ');

PlotPL('lipolysis', 0)
PlotUncertainty(0,res, 'lipolysis', 0)