% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% --------------------- Curved Microtubule Detector ----------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

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

clear; clc; close all;
addpath( genpath(pwd) )

% Define the nd2 movie that will be analyzed
filepath = '/Users/saadjansari/Documents/Projects/FY Datasets/';
filename = '998_150msR_100G_trig_7Z_001';

% Extract individual cells from an nd2 movie. The cells are saved in segmented form and can be loaded later
extractCellsFromMovie( [filepath, filename, '.nd2'] );
cellpath = [filepath, filename, '.mat'];

% Analysis will be performed by cell and by time (different cells can be run on different cores)
% for each cell and at given time, a decision will be made on the mitotic state of the cell( interphase or  prophase/metaphase). For now, this can be based on the degree of localization of intensity, i.e how well is intensity distributed inside the cell( highly concentrated or spread out). Also a bias will be added so the state of the cell is most likely whatever the state was at time t-1. Later on, a CNN-based approach could be used to classify the cell phase.

% How is all this laid out? 
% Cell Loop Begin
% Detection Step( time loop begin, cellphase find, estimate features, fit features locally, fit features globally, add/remove features, fit globally again, time loop end)
% in the detection step, all data needed for plotting is saved( say in a csv) and an isolated plotting package is called to generate the plots (so that things can be run on a remote server and plotting can be performed locally after data is pulled down.)
% Tracking Step( 

% Load posit file defining parameters( interphase, metaphase, anaphase, kc etc)
params = 1;
mov = matfile( cellpath, 'Writable', true);

for jCell = 1 : mov.NumCells

	% detect features in all time frames
	featuresRaw = detectFeaturesCell( mov, jCell, params)

	% track features over time
	featuresTracked = trackFeaturesCell( featureInfo, jCell, params)

	% run analysis on features
	featuresAnalyzed = analyzeFeaturesCell( featuresTracked, jCell, params) 
end




% For each cell, estimate and fit microtubules
for jCell = 1 : length( IsolatedCells)
    
    currCell = IsolatedCells(jCell).cellNumber; % current cell number
    savePathCell = [ savePath, filesep, 'cell_', num2str(currCell)]; % path to save any/all results in

    % Extract correct cell information
    imCellTimeline = IsolatedCells( jCell).cell3D;
    imCellTimelineRaw = IsolatedCells( jCell).raw;
    
    for jFrame = 1 : length( framesToAnalyze)
        
        currFrame = framesToAnalyze( jFrame); % current frame
        
        savePathCellFrame = [ savePathCell, filesep, 'frame_' num2str(currFrame) ];
        mkdir( savePathCellFrame); % create this path
        addpath( savePathCellFrame) % add path to matlab
        diary( [savePathCellFrame, filesep, 'diary.txt'] ) % start diary to record output progression.
    
        disp( ['Movie Name: ' metaData.fileName])
        disp( ['Current Cell: ' num2str(currCell)])
        disp( ['Current Frame: ' currFrame])
        disp( ['Data Location: ' savePathCellFrame])
        
        % frames to use for estimation. Use 3,5 frames for better results.
        frameRange = 5; % use odd number
        if currFrame <= ceil( (frameRange-1) / 2 )
            framesforGuess = 1 :  currFrame + ceil( (frameRange-1) / 2 );
        elseif currFrame > size( imCellTimeline,3) - ceil( (frameRange-1) / 2 )
            framesforGuess = currFrame - ceil( (frameRange-1) / 2 ) : size( imCellTimeline,3);
        else
            framesforGuess = currFrame - ceil( (frameRange-1) / 2 ) : currFrame + ceil( (frameRange-1) / 2 );
        end
    
        imageForGuess = IsolatedCells( jCell).raw(:, :, framesforGuess); % image for guess
        maskForCell = logical( IsolatedCells( jCell).cell3D(:, :, currFrame) ); % mask from segmentation
        
        % Generate estimates of microtubules
        MicrotubuleBank = generateGuessInterphase( imageForGuess, maskForCell, savePathCellFrame);

        % Fit microtubules
    	imageForFit = IsolatedCells( jCell).raw(:, :, currFrame) .* maskForCell; % image for fitting
        MicrotubuleBank = fitInterphaseCurves( MicrotubuleBank, imageForFit, savePathCellFrame);
        
        % save cell number, and segmentation mask with cell location, and cell
        % movie in the folder.
        cellFrameInfo.cellNumber = currCell;
        cellFrameInfo.frameNumber = currFrame;
        cellFrameInfo.frameTotal = size( imCellTimeline,3);
        cellFrameInfo.framesSeg = framesForSeg;
        cellFrameInfo.location = IsolatedCells(jCell).locations;
        cellFrameInfo.MicrotubuleBank = MicrotubuleBank;
        save( [ savePathCellFrame, filesep, 'cellFrameInfo.mat'], 'cellFrameInfo')
        clear MicrotubuleBank cellFrameInfo
        diary off
    
    end
    
end

% for jCell = 1 : length( IsolatedCells)
%     currCell = IsolatedCells(jCell).cellNumber; % current cell number
%     savePathCell = [ savePath, filesep, 'cell_', num2str(currCell)]; % path to save any/all results in
% 
%     % Extract correct cell information
%     imCellTimeline = IsolatedCells( jCell).cell3D;
%     imCellTimelineRaw = IsolatedCells( jCell).raw;
%     currFrame = 1; % current frame
%     savePathCellFrame = [ savePathCell, filesep, 'frame_' num2str(currFrame) ];
%     frameRange = 5; % use odd number
%     if currFrame <= ceil( (frameRange-1) / 2 )
%         framesforGuess = 1 :  currFrame + ceil( (frameRange-1) / 2 );
%     elseif currFrame > size( imCellTimeline,3) - ceil( (frameRange-1) / 2 )
%         framesforGuess = currFrame - ceil( (frameRange-1) / 2 ) : size( imCellTimeline,3);
%     else
%         framesforGuess = currFrame - ceil( (frameRange-1) / 2 ) : currFrame + ceil( (frameRange-1) / 2 );
%     end
%     imageForGuess = IsolatedCells( jCell).raw(:, :, framesforGuess); % image for guess
%     maskForCell = logical( IsolatedCells( jCell).cell3D(:, :, currFrame) ); % mask from segmentation
%     MicrotubuleBank = generateGuessInterphase( imageForGuess, maskForCell, savePathCellFrame);
%     clear MicrotubuleBank
% end
