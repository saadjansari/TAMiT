function [CFG, CFGinfo] = initConfiguration( CFG)
    % Configuration File: sets up the flags and the paths for running the software in different envoronments

    % Configurations: { Local, Debug, Summit, Rumor}
    if nargin == 0
        CFG = 'Local';
    end

    if ~strcmp(CFG, 'Local') &&  ~strcmp(CFG, 'Debug') && ~strcmp(CFG, 'Summit') && ~strcmp(CFG, 'Rumor')  
        error('makeSettings: CFG input must be either ''Local'', ''Debug'', ''Summit'', or ''Rumor'' ');
    end

    % Initialize flags
    runFit = 0;
    runRealTimeGraphics = 0;
    runPostFitGraphics = 0;
    runTracking = 0;
    runAnalysis = 0;
    runDebug = 0;
    fitLocal = 0;
    fitFeatureNumber = 1;
    fitExploreSpeed = 0;

    switch CFG
        case 'Local'
            runFit = 1;
            runRealTimeGraphics = 0;
            runPostFitGraphics = 0;
            runTracking = 1;
            runAnalysis = 1;
            runDebug = 0;

        case 'Debug'
            runFit = 1;
            runRealTimeGraphics = 1;
            runPostFitGraphics = 1;
            runTracking = 1;
            runAnalysis = 1;
            runDebug = 1;
            fitLocal = 1;
            fitFeatureNumber = 1;

        case 'Summit' 
            runFit = 1;
            runRealTimeGraphics = 0;
            runPostFitGraphics = 1;
            runTracking = 1;
            runAnalysis = 1;
            runDebug = 0;
            fitLocal = 1;
            fitFeatureNumber = 1;

        case 'Rumor'
            runFit = 0;
            runRealTimeGraphics = 0;
            runPostFitGraphics = 0;
            runTracking = 0;
            runAnalysis = 0;
            runDebug = 0;
            fitLocal = 1;
            fitFeatureNumber = 1;

    end

    % Forced flags
    if runDebug, runRealTimeGraphics=1; end
    if runDebug, runPostFitGraphics=1; end

    % Assign paths
    % runPath : location where the main.m file will be run from.
    % saveParent : the parent directory where the results folder will reside
    if strcmp( CFG, 'Local')

        runPath = '/Users/saadjansari/Documents/Projects/ImageAnalysis/SingleCell';
        saveParent = '/Users/saadjansari/Documents/Projects/ImageAnalysis/SingleCell/Results';

    elseif strcmp( CFG, 'Debug')

        runPath = '/Users/saadjansari/Documents/Projects/ImageAnalysis/SingleCell';
        saveParent = '/Users/saadjansari/Documents/Projects/ImageAnalysis/SingleCell/Results';

    elseif strcmp( CFG, 'Summit') 

        runPath = '/projects/saan8193/ImageAnalysis/SingleCell';
        saveParent = '/scratch/summit/saan8193/SingleCell';

    elseif strcmp( CFG, 'Rumor') 

        error('Rumor not setup for singleCell yet')
        runPath = '/projects/saan8193/ImageAnalysis/SingleCell';
        saveParent= '/scratch/summit/saan8193/ImageAnalysis/SingleCell/Results'

    end

    CFGinfo.runFit = runFit;
    CFGinfo.runRealTimeGraphics = runRealTimeGraphics;
    CFGinfo.runPostFitGraphics = runPostFitGraphics;
    CFGinfo.runTracking = runTracking; 
    CFGinfo.runAnalysis = runAnalysis; 
    CFGinfo.runDebug = runDebug; 
    CFGinfo.runPath = runPath; 
    CFGinfo.saveParent = saveParent; 
    CFGinfo.fitLocal = fitLocal;
    CFGinfo.fitFeatureNumber = fitFeatureNumber;
    CFGinfo.fitExploreSpeed = fitExploreSpeed;

end
