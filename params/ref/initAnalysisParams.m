function params = initAnalysisParams();

    % default parent path where fitting results will exist, unless a full path is provided to analysis function
    params.pathParent = '/Users/saadjansari/Documents/Projects/ImageAnalysis/SingleCell/Results';

    % Cell type. Mitotic analysis differs from interphase analysis.
    params.cellType = 'Monopolar';

    % channelsToAnalyze gives the channel index where channelFeatures are located
    params.channelsToAnalyze = [ 2, 1];
    params.channelFeatures = {'Microtubule', 'Cut7'};
    params.timeStep = 8; % Can be overwritten in analysis code

    params.flagMovie = 0;
    params.flagGraph = 0;

end
