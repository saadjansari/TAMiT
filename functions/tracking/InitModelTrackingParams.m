function TrackParams = InitModelTrackingParams( filename, params)
% InitModelTrackingParams : Initializes tracking parameters for each model
% used for feature detection

    % Defining the parameters
    % Gap Closing General (gapCloseParams)
        % timeWindow : Max time between track segment end and track segment start for successful linking
        % mergeSplit : Merging and Splitting ( 1 for both allowed, 2 for merging allowed, 3 for splitting allowed, 0 for neither allowed)
        % minTrackLen : Minimum length of track segments to be used for closing gaps (in frames)
        % diagnostics : 1 to plot a histrogram of gap lengths, 0 otherwise (optional argument)
    % Cost Matrices Linking: 
        % funcName : name of function used in cost calculation
        % parameters
            % linearMotion : Type of Motion (0 for simple diffusion only, 1 for simple diffusion + directed motion, 2 for simple diffusion + switchable directed motion)
            % minSearchRadius : min search radius (determines how far an object could have travelled betwen frames) (lower bound)
            % maxSearchRadius : max search radius (determines how far an object could have travelled betwen frames) (upper bound)
            % brownStdMult : multiplication factor (multiplied with standard deviation of positions to calculate search radius)
            % useLocalDensity : flag for using local densities to expand the search radius of isolated features.
            % nnWindow : number of frames to use to determine isolation of feature in time.
            % diagnostics : max frame-frame distance to include in histogram. Leave empty for no histogram (optional argument).
    % Cost Matrices Gap Closing: 
        % funcName : name of function used in cost calculation
        % parameters
            % linearMotion : Type of Motion (0 for simple diffusion only, 1 for simple diffusion + directed motion, 2 for ?????) 
            % minSearchRadius : min search radius (determines how far an object could have travelled betwen frames) (lower bound)
            % maxSearchRadius : max search radius (determines how far an object could have travelled betwen frames) (upper bound)
            % brownStdMult : multiplication factor (multiplied with standard deviation of positions to calculate search radius for brownian motion)
            % diagnostics : max frame-frame distance to include in histogram. Leave empty for no histogram (optional argument).
            % timeReachConfB : time after which to change multiplicative factor for search radius for brownian motion
            % brownScaling : power for scaling brownian search radius before and after timeReachConfB
            % ampRatioLimit : ampitude ratio for intensity before and after merging, or, after and before splitting. Ideally this should be 1.
            % lenForClassify : minimum track length for any type of motion. we expect that a particle undergoes a type of motion for atleast this many frames
            % useLocalDensity : flag for using local densities to expand the search radius of isolated features.
            % nnWindow : number of frames to use to determine isolation of feature in time.
            % linStdMult : multiplication factor (multiplied with standard deviation of tracks to calculate search radius for linear motion)
            % timeReachConfL : time after which to change multiplicative factor for search radius of linear directed motion
            % linScaling : power for scaling linear search radius before and after timeReachConfL
            % maxAngleVV : max difference in angle for closing two tracks
            % gapPenalty : penalty for gaps, 1 for no penalty. particle disappearance for N frames leads to penalty of gapPen^N (optimal argument)
            % resLimit : resolution limit, generally equal to 3*PSF

    % Gap-Closing Parameters {{{
    gapCloseParam.timeWindow = 5;
    gapCloseParam.mergeSplit = 1; 
    gapCloseParam.minTrackLen = 2;
    gapCloseParam.diagnostics = 1;
    % }}}

    % Cost Matrix for frame-to-frame linking {{{
    costMatrices(1).funcName = 'costMatRandomDirectedSwitchingMotionLink_Extended';
    parameters.linearMotion = 0;
    parameters.minSearchRadius = 1; 
    parameters.maxSearchRadius = 10;
    parameters.brownStdMult = 4;
    parameters.useLocalDensity = 1; 
    parameters.nnWindow = gapCloseParam.timeWindow;
    % Kalman Filter Initialization Parameters
    parameters.kalmanInitParam = [];
    parameters.kalmanInitParam.searchRadiusFirstIteration = 5;
    parameters.diagnostics = 30; 
    costMatrices(1).parameters = parameters;
    clear parameters
    % }}}

    % Cost Matrix for Gap Closing {{{
    costMatrices(2).funcName = 'costMatRandomDirectedSwitchingMotionCloseGaps';
    parameters.linearMotion = 1;
    parameters.minSearchRadius = 1; 
    parameters.maxSearchRadius = 10;
    parameters.brownStdMult = 3*ones(gapCloseParam.timeWindow,1);
    parameters.timeReachConfB = gapCloseParam.timeWindow;
    parameters.brownScaling = [0 0.01];
    parameters.ampRatioLimit = [0.25 4];
    parameters.lenForClassify = 5;
    parameters.useLocalDensity = 1; 
    parameters.nnWindow = gapCloseParam.timeWindow;
    parameters.linStdMult = 4*ones(gapCloseParam.timeWindow,1);
    parameters.timeReachConfL = gapCloseParam.timeWindow;
    parameters.linScaling = [0.25 0.01];
    parameters.maxAngleVV = 45;
    parameters.gapPenalty = 1;
    parameters.resLimit = []; 
    costMatrices(2).parameters = parameters;
    clear parameters
    % }}}

    % Other Parameters {{{
    % Kalman filter function names
    kalmanFunctions.reserveMem  = 'kalmanResMemLM';
    kalmanFunctions.initialize  = 'kalmanInitLinearMotion';
    kalmanFunctions.calcGain    = 'kalmanGainLinearMotion';
    kalmanFunctions.timeReverse = 'kalmanReverseLinearMotion';
    % additional input
    savedir = params.savePath;
    filenamebase = filename;
    %saveResults.dir = savedir; %directory where to save input and output
    %saveResults.filename = [filenamemain,'_cell_',num2str(cellNum),'_tracks.mat']; %name of file where input and output are saved
    saveResults = 0; %don't save results
    %verbose state
    verbose = 1;
    %problem dimension
    probDim = 2;
    % }}}

    TrackParams.gapCloseParam = gapCloseParam;
    TrackParams.costMatrices = costMatrices;
    TrackParams.kalmanFunctions = kalmanFunctions;
    TrackParams.saveResults = saveResults;
    TrackParams.verbose = verbose;
    TrackParams.probDim = probDim;

end
