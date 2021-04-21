function params = initAnalysisParams();

    % default parent path where fitting results will exist, unless a full path is provided to analysis function
%     params.pathParent = '/Users/saadjansari/Documents/Projects/ImageAnalysis/SingleCell/Results/Test3';
%     params.pathParent = '/Users/saadjansari/Documents/Projects/ImageAnalysis/SingleCell/Results/Paper/Bipolar';
    params.pathParent = '/Users/saadjansari/Documents/Projects/ImageAnalysis/SingleCell/Results/Paper/Monopolar';
%         params.pathParent = '/Users/saadjansari/Documents/Projects/ImageAnalysis/SingleCell/Results/Paper/fig3/pombe_anaphase_elongation';

    % Cell type. Mitotic analysis differs from interphase analysis.
    params.cellType = 'MitosisBud';

    % channelsToAnalyze gives the channel index where channelFeatures are located
    params.channelsToAnalyze = [ 1];
    params.channelFeatures = {'Microtubule'};
    params.timeStep = 8; % Can be overwritten in analysis code

    params.flagMovie = 1;
    params.flagGraph = 0;

end
