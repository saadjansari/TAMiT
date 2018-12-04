function featureInfo = detectFeaturesCell( movData, currentCell, params)
% explanation needed

msg = { ['Analyzing Cell-', num2str(currentCell)] };disp([msg{:}])

% Load metadata
params.meta = movData.metaData;
params.movieMatfile = movData;
savePath = params.savePath;
% cellphase vector will store cellphase info (1 is interphase, 2 is mitosis, 3 is post-anaphase)
cellphase = zeros( 1, params.meta.numTimes);

times = 1:300;
channels = 1;

% Pre-allocate featureBank {{{
paramsInit = params; paramsInit.savePath = pwd; paramsInit.currTime = 1; paramsInit.currChannel = 1;
ftB = initFeatureBank( 0,0, 1, 1, 1, paramsInit);
for jT = 1 : times(end), for jC = 1 : channels(end) 
    featureInfo(jT, jC) = ftB;
end; end
size( featureInfo)
% }}}

for jTime = times 
for jChannel = channels 
  
    params.savePath = [savePath, filesep, sprintf('cell_%d',currentCell), filesep, sprintf('channel_%d', jChannel), filesep, sprintf('frame_%d', jTime) ]; 
    params.currTime  = jTime;
    params.currChannel = jChannel;
    mkdir( params.savePath)

	% load image for estimation
	[ imForEstimate, imForEstimate_mt, mask] = loadImageForEstimation( movData, jTime, jChannel, params.estimation.InitImageFrameRange);

	% classify cell phase
	cellphase = classifyCellPhase( imForEstimate.*mask, cellphase, params.cellPhaseClassification);

	% Initialize feature bank based on cellphase
	featureInfo( jTime, jChannel) = initFeatureBank( imForEstimate, mask, currentCell, cellphase(jTime), mean( movData.planeTimes( :, jTime, jChannel) ), params);

	% Estimate features
 	featureInfo( jTime, jChannel) = estimateFeatures( featureInfo, imForEstimate_mt, cellphase(jTime), params);

	% Fit features
  	featureInfo( jTime, jChannel) = fitFeatures( featureInfo( jTime, jChannel), cellphase(jTime), params);
    
end
end
    
save( [savePath, filesep, 'featureInfo.mat'], 'featureInfo', 'params');

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
function ftBank = initFeatureBank( imPlane, mask, currentCell, currentCellPhase, currentTime, initParams)
	

	switch currentCellPhase
	case 1
        initParams = rmfield( initParams, {'mitosis', 'kc'} );
		ftBank = feval( initParams.interphase.func.initializeBank, imPlane, mask, currentCell, currentTime, initParams);
	case 2
        initParams = rmfield( initParams, {'interphase', 'kc'} );
		ftBank = feval( initParams.mitosis.func.initializeBank, imPlane, initParams);
	otherwise
        initParams = rmfield( initParams, {'mitosis', 'kc'} );
		ftBank = feval( initParams.interphase.func.initializeBank, imPlane, initParams);
    end
    

end
%j }}}

% obj = estimateFeatures( imPlane, currentCellPhase, estParams) {{{
function featureBank = estimateFeatures( featureInfo, imXYZT, currentCellPhase, estParams)

    % check if estimation process should be run again or if old fitted data should be taken as the estimate for the next time frame
    if estParams.estimation.usePreviousFitAsEstimate && params.currTime  > times(1) 
        useFitAsEstimate = 1; else, useFitAsEstimate = 0; end

	switch currentCellPhase
	case 1 

    % Interphase {{{
        estParams = rmfield( estParams, {'mitosis', 'kc'} );
        if useFitAsEstimate
            oldTime = times( find( times == params.currTime) - 1);
		    featureBank = feval( estParams.interphase.func.estimateFeaturesFromPrior, featureInfo( params.currTime, params.currChannel), featureInfo( oldTime, params.currChannel), imXYZT, estParams);
        else
		    featureBank = feval( estParams.interphase.func.estimateFeaturesDeNovo, featureInfo( params.currTime, params.currChannel), imXYZT, estParams);
        end
    %  }}}

	case 2

    % Mitosis {{{
        estParams = rmfield( estParams, {'interphase', 'kc'} );
        if useFitAsEstimate
            oldTime = times( find( times == params.currTime) - 1);
		    featureBank = feval( estParams.interphase.func.estimateFeaturesFromPrior, featureInfo( params.currTime, params.currChannel), featureInfo( oldTime, params.currChannel), imXYZT, estParams);
        else
		    featureBank = feval( estParams.interphase.func.estimateFeaturesDeNovo, featureInfo( params.currTime, params.currChannel), imXYZT, estParams);
        end
		featureBank = feval( estParams.mitosis.func.estimateFeatures, featureBank, estParams);
    %  }}}

	otherwise

    %  Otherwise {{{
        estParams = rmfield( estParams, {'mitosis', 'kc'} );
        if useFitAsEstimate
            oldTime = times( find( times == params.currTime) - 1);
		    featureBank = feval( estParams.interphase.func.estimateFeaturesFromPrior, featureInfo( params.currTime, params.currChannel), featureInfo( oldTime, params.currChannel), imXYZT, estParams);
        else
		    featureBank = feval( estParams.interphase.func.estimateFeaturesDeNovo, featureInfo( params.currTime, params.currChannel), imXYZT, estParams);
        end
		featureBank = feval( estParams.interphase.func.estimateFeatures, featureBank, estParams);
	end
    % }}}

end
 % }}}

% obj = fitFeatures( imPlane, currentCellPhase, fitParams) {{{
function featureBank = fitFeatures( featureBank, currentCellPhase, fitParams)

	switch currentCellPhase
	case 1 
        fitParams.interphase.savePath = fitParams.savePath;
        fitParams.interphase.config = fitParams.config;
		featureBank = feval( fitParams.interphase.func.fitFeatures, featureBank, fitParams.interphase);
	case 2
        fitParams.mitosis.savePath = fitParams.savePath;
        fitParams.mitosis.config = fitParams.config;
		featureBank = feval( fitParams.mitosis.func.fitFeatures, featureBank, fitParams.mitosis);
	otherwise
        fitParams.interphase.savePath = fitParams.savePath;
        fitParams.interphase.config = fitParams.config;
		featureBank = feval( fitParams.interphase.func.fitFeatures, featureBank, fitParams.interphase);
	end

end
 % }}}


end
