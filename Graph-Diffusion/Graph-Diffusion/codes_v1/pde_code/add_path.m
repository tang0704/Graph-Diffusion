function add_path
% Add folders of predefined functions into matlab searching paths.

global footpath;
footpath = cd;

addpath(genpath([footpath '/lib']));
addpath(genpath([footpath '/img']));
addpath(genpath([footpath '/alg']));
addpath(genpath([footpath '/evaluate']));

addpath(genpath([footpath '/Code_LowRankSaliency']));