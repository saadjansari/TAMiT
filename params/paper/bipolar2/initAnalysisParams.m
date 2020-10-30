function params = initAnalysisParams();

    % default parent path where fitting results will exist, unless a full path is provided to analysis function
    params.pathParent = '/Users/saadjansari/Documents/Projects/ImageAnalysis/SingleCell/Results/Paper';

    % Cell type. Mitotic analysis differs from interphase analysis.
    params.cellType = 'Mitosis';

    % channelsToAnalyze gives the channel index where channelFeatures are located
    params.channelsToAnalyze = [ 1];
    params.channelFeatures = {'Microtubule'};
    params.timeStep = 8; % Can be overwritten in analysis code

    params.flagMovie = 1;
    params.flagGraph = 0;

end
