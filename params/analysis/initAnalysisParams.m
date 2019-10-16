function params = initAnalysisParams();

    % default parent path where fitting results will exist, unless a full path is provided to analysis function
    params.pathParent = '/Users/saadjansari/Documents/Projects/ImageAnalysis/SingleCell/Results';

    % Cell type. Mitotic analysis differs from interphase analysis.
    params.cellType = 'Mitosis';

    % channelsToAnalyze gives the channel index where channelFeatures are located
    params.channelsToAnalyze = [ 1, 2];
    params.channelFeatures = {'Cut7', 'Microtubule'};
    params.timeStep = 8; % Can be overwritten in analysis code

    params.flagMovie = 1;
    params.flagGraph = 1;

end
