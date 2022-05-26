function params = initAnalysisParams()

    % default parent path where fitting results will exist, unless a full path is provided to analysis function
    params.pathParent = '/Users/saadjansari/Documents/Projects/ImageAnalysis/SingleCell/Results';
    fprintf('NOTE: User must set correct path for Analysis!\n')
    fprintf('Current Path: %s', params.pathParent)

    % channelsToAnalyze gives the channel index where channelFeatures are located
    params.channelsToAnalyze = [ 1,2];
    params.channelFeatures = {'Microtubule','Cut7'};
    params.timeStep = 8; % Can be overwritten in analysis code

    params.flags.movie = 0;
    params.flags.graph = 0;
    params.flags.tracking = 1;
    params.flags.tracking_movie = 1;
    
    % Paper Figure flags
    params.flags.paper_figs = 0;

end