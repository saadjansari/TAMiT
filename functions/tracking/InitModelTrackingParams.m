function TrackParams = InitModelTrackingParams( filename)
% InitModelTrackingParams : Initializes tracking parameters for each model
% used for feature detection

    %% General Gap-Closing Parameters
    
    % Max time between track segment end and track segment start for
    % successful linking
    gapCloseParam.timeWindow = 5;
    
    % Merging and Splitting ( 1 for both allowed, 2 for merging allowed, 3
    % for splitting allowed, 0 for neither allowed)
    gapCloseParam.mergeSplit = 1; 
    
    % Minimum length of track segments to be used for closing gaps (in
    % frames)
    gapCloseParam.minTrackLen = 2;

    % 1 to plot a histrogram of gap lengths, 0 otherwise (optional
    % argument)
    gapCloseParam.diagnostics = 1;
    
    %% Cost Matrix for frame-to-frame linking
    
    % Function Name
    costMatrices(1).funcName = 'costMatRandomDirectedSwitchingMotionLink';
    
    % Parameters
    % Type of Motion : 0 for simple diffusion only, 1 for simple diffusion
    % + directed motion 
    parameters.linearMotion = 2;
    
    % Min Search Radius ( The min radius that determines how far an object
    % could have travelled betwen frames) (calculated in code but bounded
    % from below by this value)
    parameters.minSearchRadius = 1; 
    
    % Max Search Radius ( The max radius that determines how far an object
    % could have travelled betwen frames) (calculated in code but bounded
    % from above by this value)
    parameters.maxSearchRadius = 10;
    
    % Multiplication factor : multiplied with standard deviation to get
    % search radius
    parameters.brownStdMult = 4;
    
    % Flag for using local densities to expand the search radius of
    % isolated features.
    parameters.useLocalDensity = 1; 

    % Number of frames to use to determine isolation of feature in time.
    parameters.nnWindow = gapCloseParam.timeWindow;
    
    % Kalman Filter Initialization Parameters
    parameters.kalmanInitParam = [];
    parameters.kalmanInitParam.searchRadiusFirstIteration = 5;

    % Linking Distances( value signifies the max frame-frame distance to
    % include in histogram). Leave empty for no histogram. Does not work
    % for the first or last frame of a movie ( optional argument).
     parameters.diagnostics = 30; 
    
    costMatrices(1).parameters = parameters;
    clear parameters
    
    %% Cost Matrix for Gap Closing

    % Function Name
    costMatrices(2).funcName = 'costMatRandomDirectedSwitchingMotionCloseGaps';
    
    % Parameters
    % Type of Motion : 0 for simple diffusion only, 1 for simple diffusion
    % + directed motion 
    parameters.linearMotion = 1;
    
    % Min Search Radius ( The min radius that determines how far an object
    % could have travelled betwen frames) (calculated in code but bounded
    % from below by this value)
    parameters.minSearchRadius = 1; 
    
    % Max Search Radius ( The max radius that determines how far an object
    % could have travelled betwen frames) (calculated in code but bounded
    % from above by this value)
    parameters.maxSearchRadius = 10;
    
    % Multiplication factor : multiplied with standard deviation to get
    % search radius. Each entry in  vector signifies the multpicative
    % constant for distance of frames 1 - gapCloseParam.timeWindow
    parameters.brownStdMult = 3*ones(gapCloseParam.timeWindow,1);
    
    % Time after which to change the stdev multiplicative factor for search
    % radius for browninan dynamics
    parameters.timeReachConfB = gapCloseParam.timeWindow;
    
    % Power for scaling the brownian search radius before and after
    % timeReachConfB
    parameters.brownScaling = [0 0.01];
    
    % Amplitude ratio for particle intesities before merging/ after
    % merging, or before splitting/after splitting. Expect this to be 1 but
    % specify a minimum and a maximum value here
    parameters.ampRatioLimit = [0.25 4];
    
    % Minimum track length for any type of motion. We will expect that a
    % particle undergoes a type of motion for atleast this many frames.
    parameters.lenForClassify = 5;
    
    % Flag for using local densities to expand the search radius of
    % isolated features.
    parameters.useLocalDensity = 1; 
    
    % Number of frames to use to determine isolation of feature in time.
    parameters.nnWindow = gapCloseParam.timeWindow;
    
    % Multiplicative factor for determining scaling of search radius for
    % linear model using standard deviation of tracks.
    parameters.linStdMult = 4*ones(gapCloseParam.timeWindow,1);
    
    % Time after which to change the stdev multiplicative factor for search
    % radius for lienar directed motion
    parameters.timeReachConfL = gapCloseParam.timeWindow;
    
    % Power for scaling the linear search radius before and after
    % timeReachConfL
    parameters.linScaling = [0.25 0.01];
    
    % Max Difference in angle of linear motion between two tracks that may
    % be closed.
    parameters.maxAngleVV = 45;
    
    % Penalty for gaps. 1 for no penalty. Particle disappearance for N
    % frames leads to a penalty of (gapPenalty)^N (optional argument)
    parameters.gapPenalty = 1;
    
    %optional; to calculate MS search radius
    %if not input, MS search radius will be the same as gap closing search radius
    parameters.resLimit = []; %resolution limit, which is generally equal to 3 * point spread function sigma.

    costMatrices(2).parameters = parameters;
    clear parameters

    %% Other Parameters

    % Kalman filter function names
    kalmanFunctions.reserveMem  = 'kalmanResMemLM';
    kalmanFunctions.initialize  = 'kalmanInitLinearMotion';
    kalmanFunctions.calcGain    = 'kalmanGainLinearMotion';
    kalmanFunctions.timeReverse = 'kalmanReverseLinearMotion';

    % additional input
    %saveResults

    savedir = [pwd '\'];
    filenamebase = filename;
    
    %saveResults.dir = savedir; %directory where to save input and output
    %saveResults.filename = [filenamemain,'_cell_',num2str(cellNum),'_tracks.mat']; %name of file where input and output are saved
    saveResults = 0; %don't save results

    %verbose state
    verbose = 1;
    %problem dimension
    probDim = 3;

    TrackParams.gapCloseParam = gapCloseParam;
    TrackParams.costMatrices = costMatrices;
    TrackParams.kalmanFunctions = kalmanFunctions;
    TrackParams.saveResults = saveResults;
    TrackParams.verbose = verbose;
    TrackParams.probDim = probDim;




end