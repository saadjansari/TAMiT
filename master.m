% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% --------------------- Curved Microtubule Detector ----------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

clear; clc; close all;
addpath( genpath(pwd) )

% run config file (runpath, filepath, savepath)
config = config_file();

% Description:
% This script can be used to segment a 2D field of view .nd2 movie of fission
% yeast cells with mCherry-atb2 tubulin and segment them. For one or more
% of those segmented cells, it can run a detection routine to detect curved
% microtubules. The detection routine begins with an estimation phase,
% where MTOCs (microtubule organizing centers) are identified as regions of
% maximum brightness inside cells. Estimates for microtubules are generated
% by starting at the MTOC and propagating at the optimal angle in a
% radially integrated map. Each MTOC can sprout 2 microtubules in this
% phase. Once the estimates are generated, a quadratic or a cubic
% polynomial is fitted through the points. The polynomial is then
% numerically convolved with a gaussian and parameters are optimized via
% 'lsqnonlin' or 'particleswarm' with respect to the original image.
% This is a standalone script implemented in MATLAB.   
%
% 1) Load the Image
%
% 2) Segmentation
%   2.1) Generate segmentation mask
%   2.2) Use segmentation mask to isolate cells
%
% 3) Estimate location of microtubules
%
% 4) Fit gaussian microtubules
%
% 5) Save any relevant information

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

% Define the nd2 movie that will be analyzed
filename = '998_150msR_100G_trig_7Z_001';

savePath = createSaveDirectory( config.savepath, filename, 1);
cellpath = [config.filepath, filename, '.mat'];

% Extract individual cells from an nd2 movie. The cells are saved in segmented form and can be loaded later
extractCellsFromMovie( [config.filepath, filename, '.nd2'] , 'all');

% Analysis will be performed by cell and by time (different cells can be run on different cores)
% for each cell and at given time, a decision will be made on the mitotic state of the cell( interphase or  prophase/metaphase). For now, this can be based on the degree of localization of intensity, i.e how well is intensity distributed inside the cell( highly concentrated or spread out). Also a bias will be added so the state of the cell is most likely whatever the state was at time t-1. Later on, a CNN-based approach could be used to classify the cell phase.

% How is all this laid out? 
% Cell Loop Begin
% Detection Step( time loop begin, cellphase find, estimate features, fit features locally, fit features globally, add/remove features, fit globally again, time loop end)
% in the detection step, all data needed for plotting is saved( say in a csv) and an isolated plotting package is called to generate the plots (so that things can be run on a remote server and plotting can be performed locally after data is pulled down.)
% Tracking Step( 

% Load posit file defining parameters( interphase, metaphase, anaphase, kc etc)
params = paramsInitialize(); % script for defining non-tunable parameters
params.savePath = savePath;
params.config = config;

mov = matfile( cellpath, 'Writable', true);

goodCells = sort([10, 11, 12, 15, 16 18, 19, 23, 24, 25, 27, 29, 31, 33, 35, 36, 37, 39, 4, 5, 7, 8]);
for jCell = 31 

	% detect features in all time frames
	featuresCell = detectFeaturesCell( mov, jCell, params);

	% track features over time
%	featuresTracked = trackFeaturesCell( featureInfo, jCell, params)

	% run analysis on features
%	featuresAnalyzed = analyzeFeaturesCell( featuresTracked, jCell, params) 
end


