function[options] = load_libraries()

addpath('~/Documents/MATLAB/minFunc_2012')
addpath('~/Documents/MATLAB/minFunc_2012')
addpath('~/Documents/MATLAB/minFunc_2012/minFunc')
addpath('~/Documents/MATLAB/minFunc_2012/minFunc/compiled')
addpath('~/Documents/MATLAB/minFunc_2012/minFunc/mex')

maxFunEvals = 10000;
options.MaxIter = 10000;
options = [];
options.GradObj = 'on';
options.Display = 'off'; %'off';%'iter';
options.useMex = 0;
options.maxFunEvals = maxFunEvals;
options.Method = 'lbfgs';