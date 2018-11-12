function featureInfo = detectFeaturesCell( movData, currentCell, params)
% explanation needed

msg = { ['Analyzing Cell-', num2str(currentCell)] };disp([msg{:}])

% Load metadata
params.meta = movData.metaData;
params.movieMatfile = movData;

% cellphase vector will store cellphase info (1 is interphase, 2 is mitosis, 3 is post-anaphase)
cellphase = zeros( 1, params.meta.numTimes);


for jTime = 1 : 1 
for jChannel = 1 : 1 
    
    if ~exist( 'tempSave.mat')
	% load image for estimation
	[ imForEstimate, imForEstimate_mt, mask] = loadImageForEstimation( movData, jTime, jChannel, params.estimation.InitImageFrameRange);

	% classify cell phase
	cellphase = classifyCellPhase( imForEstimate.*mask, cellphase, params.cellPhaseClassification);

	% Initialize feature bank based on cellphase
	featureBank = initFeatureBank( imForEstimate, mask, currentCell, cellphase(jTime), mean( movData.planeTimes( :, jTime, jChannel) ) ,params);

	% Estimate features
 	featureBank = estimateFeatures( featureBank, imForEstimate_mt, cellphase(jTime), params);
    save( 'tempSave.mat', 'featureBank', 'cellphase'); 
    else
    load('tempSave.mat')
    end
	% Fit features
  	featureBank = fitFeatures( featureBank, cellphase(jTime), params)
    save('tempSaveFit.mat', 'featureBank', 'cellphase', 'params');

    featureInfo( jTime, jChannel) = featureBank;

end
end

% CellPhasePlot( movData, cellphase);


% Functions

% imEstXYZ = loadImageForEstimation( movMatFile,currentTime,currentChannel,frameRange) {{{
function [imEstXYZ, imEstXYZT, mask] = loadImageForEstimation( movMatFile,currentTime,currentChannel,frameRange)

	if currentTime  <= ceil( (frameRange-1) / 2 )
            framesforGuess = 1 :  currentTime + ceil( (frameRange-1) / 2 );
        elseif currentTime  > params.meta.numTimes - ceil( (frameRange-1) / 2 )
            framesforGuess = currentTime - ceil( (frameRange-1) / 2 ) : meta.numTimes;
        else
            framesforGuess = currentTime - ceil( (frameRange-1) / 2 ) : currentTime + ceil( (frameRange-1) / 2 );
        end

	imEstXYZT = mat2gray( movMatFile.(['cellRaw_', num2str(currentCell)] )(:,:,:,framesforGuess, currentChannel) );
	imEstXYZ = squeeze( mean( imEstXYZT, 4) );
    mask = logical( movMatFile.(['cell_', num2str(currentCell)] )(:,:,:,framesforGuess, currentChannel) );
    mask = logical( mean( mask, 4) );

end
% }}}

% cellphase = classifyCellPhase( imPlane, cellphase, classificationParams) {{{
function cellphase = classifyCellPhase( imPlane, cellphase, classificationParams)

	% will be based on the tubulin channel
	if jChannel == 2 % only do this for tubulin channel
		return
	end

	% we care about the distribution of intensity across pixels inside the cell
	imPlane = max( imPlane, [], 3); % MIP in z-dimension
	imPlane( imPlane == 0) = median( imPlane( imPlane ~= 0) );
	% find pixels for which SnR is greater than some number
	pixVal = imPlane(:);
	noise = median( pixVal );
	snr = pixVal/noise; 
	highSNR = find( snr > classificationParams.mitosisSNR);	
	if length( highSNR) > classificationParams.mitosisMinPixels
		cellphase( jTime) = 2; % prophase/metaphase/anaphase
	else
		cellphase( jTime) = 1; % interphase
	end
    
    warning('detectFeaturesCell: CellPhase has been set to interphase by admin!')
    cellphase( jTime) = 1;


end 
% }}}

% CellPhasePlot( movMatFile, cellphase) {{{
function CellPhasePlot( movMatFile, cellphase)

	figure; set( gcf, 'WindowState', 'maximized');
	for jTime = 1 : meta.numTimes
		
		imPlane = max( loadImageForEstimation( movMatFile, jTime, 1, params.Estimation.InitImageFrameRange), [], 3);
		if cellphase(jTime) == 2, jcellphase = 'mitosis'; elseif cellphase(jTime) == 1, jcellphase = 'interphase'; else, jcellphase = 'TBD'; end
		
		imagesc(imPlane); colormap gray;
		title( [ 'Time ', num2str(jTime), ' : ', jcellphase] )
		drawnow
		pause(0.1)

	end


end
% }}}

% featureBank = initFeatureBank( imPlane, currentCellPhase, initParams) {{{
function featureBank = initFeatureBank( imPlane, mask, currentCell, currentCellPhase, currentTime, initParams)
	

	switch currentCellPhase
	case 1
        initParams = rmfield( initParams, {'mitosis', 'kc'} );
		featureBank = feval( initParams.interphase.func.initializeBank, imPlane, mask, currentCell, currentTime, initParams);
	case 2
        initParams = rmfield( initParams, {'interphase', 'kc'} );
		featureBank = feval( initParams.mitosis.func.initializeBank, imPlane, initParams);
	otherwise
        initParams = rmfield( initParams, {'mitosis', 'kc'} );
		featureBank = feval( initParams.interphase.func.initializeBank, imPlane, initParams);
	end


end
%j }}}

% obj = estimateFeatures( imPlane, currentCellPhase, estParams) {{{
function featureBank = estimateFeatures( featureBank, imXYZT, currentCellPhase, estParams)

	switch currentCellPhase
	case 1 
        estParams = rmfield( estParams, {'mitosis', 'kc'} );
		featureBank = feval( estParams.interphase.func.estimateFeatures, featureBank, imXYZT, estParams);
	case 2
        estParams = rmfield( estParams, {'interphase', 'kc'} );
		featureBank = feval( estParams.mitosis.func.estimateFeatures, featureBank, estParams);
	otherwise
        estParams = rmfield( estParams, {'mitosis', 'kc'} );
		featureBank = feval( estParams.interphase.func.estimateFeatures, featureBank, estParams);
	end

end
 % }}}

% obj = fitFeatures( imPlane, currentCellPhase, fitParams) {{{
function featureBank = fitFeatures( featureBank, currentCellPhase, fitParams)

	switch currentCellPhase
	case 1 
		featureBank = feval( fitParams.interphase.func.fitFeatures, featureBank, fitParams.interphase);
	case 2
		featureBank = feval( fitParams.mitosis.func.fitFeatures, featureBank, fitParams.mitosis);
	otherwise
		featureBank = feval( fitParams.interphase.func.fitFeatures, featureBank, fitParams.interphase);
	end

end
 % }}}


end
