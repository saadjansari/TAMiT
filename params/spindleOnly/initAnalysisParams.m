function params = initAnalysisParams();

    % default parent path where fitting results will exist, unless a full path is provided to analysis function
    params.pathParent = '/Users/saadjansari/Documents/Projects/ImageAnalysis/SingleCell/Results';

    % Is it imported from Summit
    params.importedFromSummit = 1; % temporary
    params.rmStringSummit = '/projects/saan8193/ImageAnalysis/FY Datasets/';
    params.addString = '/Users/saadjansari/Documents/Projects/ImageAnalysis/FY Datasets/';

    % Cell type. Mitotic analysis differs from interphase analysis.
    params.cellType = 'Mitosis';

    % channelsToAnalyze gives the channel index where channelFeatures are located
    params.channelsToAnalyze = [ 2, 1];
    params.channelFeatures = {'Microtubule', 'Cut7'};
    params.timeStep = 8; % Can be overwritten in analysis code


    % temporarily save params in the parent path
%     save( params.pathParent, params);

end
