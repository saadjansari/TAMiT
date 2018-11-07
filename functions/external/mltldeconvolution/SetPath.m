function SetPath()
% SetPath()
%
% Update the search path with all the directories
% necessary for running the MLTL algorithm.
%
% (c) Cedric Vonesch, Biomedical Imaging Group, EPFL, 2007.02.21-2009.04.16

basepath = cd();

addpath([basepath, '/Misc']);
addpath([basepath, '/MLTL']);
addpath([basepath, '/Wavelets']);
addpath([basepath, '/Wavelets/Compiled']);