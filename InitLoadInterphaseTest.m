% Script to segment, select, and detect interphase curved microtubules:
clear; clc; close all;
addpath( genpath(pwd) )

%% Image Load
reboot = 0; 

if reboot
    interphase_loadfile;
    imTub = imData(:,:,:,:,1);
    imTub = squeeze( max(imTub, 3) );
    save([pwd, filesep, 'imData'], 'imData'); 
    save([pwd, filesep, 'imTub'], 'imTub'); 
    clearvars -except imTub
else
    load('imTub.mat')
end

movieName = 'anon';

folderTimeStamp = datestr( now, 'yymmdd_HHMM');
SavePath = [ pwd, filesep, 'Pictures', filesep, folderTimeStamp];

SavePath = [ SavePath, filesep, movieName];


%% Process field of view movie (doesnt work well)

% segment cells based on a range of frames.
framesForSeg = 1:75;
imForSeg = mat2gray( sum( imTub( :, :, framesForSeg), 3 ) );
SegmentationInfo = generateSegmentationMask( imForSeg);

IsolatedCells = useSegmentationMask( imTub( :, :, framesForSeg), SegmentationInfo.MaskLogical, [1:SegmentationInfo.NumCells]);
% Cells have been isolated.

for jCell = 36 : length( IsolatedCells)
    
    currCell = IsolatedCells(jCell).cellNumber;
    SavePathCell = [ SavePath, filesep, num2str(currCell)];
    %create this path
    mkdir(SavePathCell)
    addpath(SavePathCell)
    
    diary( [SavePathCell, filesep, 'diary.txt'] )
    disp(['Movie Name: ' movieName])
    disp(['Current Cell: ' num2str(currCell)])
    disp(['Data Location: ' SavePathCell])

    %% Generate Estimates

    % Generate a guess for the microtubule locations.
    currFrame = 2;
    framesforGuess = currFrame-1 : currFrame+1;
    imCellTimeline = IsolatedCells( jCell).cell3D;
    imCellTimelineRaw = IsolatedCells( jCell).raw;

    % save cell number, and segmentation mask with cell location, and cell
    % movie in the folder.
    cellInfo.cellNumber = currCell;
    cellInfo.framesSeg = framesForSeg;
    cellInfo.location = IsolatedCells(jCell).locations;

    logicMask = logical( IsolatedCells( jCell).cell3D(:, :, currFrame) );
    imageFull = IsolatedCells( jCell).raw(:, :, framesforGuess);

    MicrotubuleBank = generateGuessInterphase( imageFull, logicMask, SavePathCell);

    %% Fit Microtubule Curves

    image2D = IsolatedCells( jCell).raw(:, :, currFrame) .* logicMask;
    MicrotubuleBank = fitInterphaseCurves( MicrotubuleBank, image2D, SavePathCell);
    
    cellInfo.MicrotubuleBank = MicrotubuleBank;
    save([SavePathCell, filesep, 'cellInfo.mat'], 'cellInfo')
    clear MicrotubuleBank cellInfo
    diary off
    
end