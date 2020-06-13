function addPath
% Add folders of predefined functions into matlab searching paths.
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 03-20-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 05-07-2013

global footpath;
footpath = cd;

addpath(genpath([footpath '/src']));
addpath(genpath([footpath '/lib']));
addpath(genpath([footpath '/conGraphs']));
addpath(genpath([footpath '/Methods']));
addpath(genpath([footpath '/utils_FM']));

addpath('./functions/')
%addpath('./svm_struct_mex/');

run('../vlfeat-0.9.20/toolbox/vl_setup');

