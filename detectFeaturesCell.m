function featureInfo = detectFeaturesCell( movData, currentCell, params)
% explanation needed

msg = { ['Analyzing Cell-', num2str(currentCell)] };disp([msg{:}])

% Load metadata
meta = movData.metaData;

% cellphase vector will store cellphase info (1 is interphase, 2 is mitosis, 3 is post-anaphase)
cellphase = zeros( 1, meta.numTimes);

featureBank = 0;

clear params;
params.Estimation.InitImageFrameRange = 5;
params.CellPhaseClassification.mitosisSNR = 3;
params.CellPhaseClassification.mitosisMinPixels = 5;

for jTime = 1 :  meta.numTimes
for jChannel = 1 : meta.numChannels

	% load image for estimation
	imForEstimate = loadImageForEstimation( movData, jTime, jChannel, params.Estimation.InitImageFrameRange);
	
	% classify cell phase
	cellphase = classifyCellPhase( imForEstimate, cellphase, params.CellPhaseClassification);

	% Initialize feature bank based on cellphase
	featureBank = initFeatureBank( imForEstimate, cellphase(jTime), params.Estimation)

	% Estimate features
	featureBank = estimateFeatures( imForEstimate, cellphase(jTime), params.Estimation)

	% Fit features
	featureBank = fitFeatures( movData, cellphase(jTime), params.Fitting)


end
end

 CellPhasePlot( movData, cellphase);

featureInfo = 1;







% Functions

% imEstXYZ = loadImageForEstimation( movMatFile,currentTime,currentChannel,frameRange) {{{
function imEstXYZ = loadImageForEstimation( movMatFile,currentTime,currentChannel,frameRange)

	if currentTime  <= ceil( (frameRange-1) / 2 )
            framesforGuess = 1 :  currentTime + ceil( (frameRange-1) / 2 );
        elseif currentTime  > meta.numTimes - ceil( (frameRange-1) / 2 )
            framesforGuess = currentTime - ceil( (frameRange-1) / 2 ) : meta.numTimes;
        else
            framesforGuess = currentTime - ceil( (frameRange-1) / 2 ) : currentTime + ceil( (frameRange-1) / 2 );
        end

	imEstXYZ = mat2gray( movMatFile.(['cell_', num2str(currentCell)] )(:,:,:,framesforGuess, currentChannel) );
	imEstXYZ = squeeze( mean( imEstXYZ, 4) );

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
function featureBank = initFeatureBank( imPlane, currentCellPhase, initParams)
	

	switch currentCellPhase
	case 1 
		featureBank = feval( initParams.interphase.func.initialize, imPlane, initParams.interphase);
	case 2
		featureBank = feval( initParams.mitosis.func.initialize, imPlane, initParams.mitosis);
	otherwise
		featureBank = feval( initParams.interphase.func.initialize, imPlane, initParams.interphase);
	end


end
% }} }

% obj = estimateFeatures( imPlane, currentCellPhase, estParams) {{{
function obj = estimateFeatures( imPlane, currentCellPhase, estParams)

	switch currentCellPhase
	case 1 
		featureBank = feval( estParams.interphase.func.estimate, imPlane, estParams.interphase);
	case 2
		featureBank = feval( estParams.mitosis.func.estimate, imPlane, estParams.mitosis);
	otherwise
		featureBank = feval( estParams.interphase.func.estimate, imPlane, estParams.interphase);
	end

end
 % }}}


% obj = fitFeatures( imPlane, currentCellPhase, estParams) {{{
function obj = fitFeatures( imPlane, currentCellPhase, fitParams)

	switch currentCellPhase
	case 1 
		featureBank = feval( estParams.interphase.func.estimate, imPlane, fitParams.interphase);
	case 2
		featureBank = feval( estParams.mitosis.func.estimate, imPlane, fitParams.mitosis);
	otherwise
		featureBank = feval( estParams.interphase.func.estimate, imPlane, fitParams.interphase);
	end

end
 % }}}
end

