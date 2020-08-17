function params = initConfiguration( opts)
    % Configuration File: sets up the flags and the paths for running the software in different envoronments

    % reverse time
    params.timeReversal = 0;
    params.newEstimateEveryT = 1;

    runFit = 1;
    runTracking = 1;
    runAnalysis = 1;
    
    switch opts.CFG
        case 'RELEASE'
            runLocalFit = 0;
            runGlobalFit = 1;
            runFeatureNumberFit = 1;
            display = 0;

        case 'DEBUG'
            runLocalFit = 0;
            runGlobalFit = 1;
            runFeatureNumberFit = 1;
            display = 1;
    end

    % Estimation {{{
    % Estimate Parameters
    estimate.display = display;

    % Spindle
    estimate.spindle.display = display;
    estimate.spindle.spindleMT = 1;
    estimate.spindle.spindlePoles = 1;
    estimate.spindle.astralMT = 1;

    % Monopolar Spindle
    estimate.monopolar.display = display;
    estimate.monopolar.pole = 1;
    estimate.monopolar.mt = 1;

    % Interphase Microtubules
    estimate.interphase.display = display;
    estimate.interphase.PolyOrder = [3 3 1];
    estimate.interphase.Visibility = 10;
    estimate.interphase.FieldOfView = 60;
    estimate.interphase.StepSize = 5;
    estimate.interphase.MinLength = 10;
    estimate.interphase.MaxCurves = 5;
    estimate.interphase.CostAcceptLink = 10;

    % Kinetochores
    estimate.kcbank.display = display;

    % Cut7
    estimate.cut7dist.display = display;
    % }}}

    % Fit {{{
    % Fit Parameters
    fit.runFit = runFit;
    fit.runLocalOptimization = runLocalFit;
    fit.runGlobalOptimization = runGlobalFit;
    fit.runFeatureNumberOptimization = runFeatureNumberFit;
    fit.useParallel = true; 
    fit.state = opts.CFG;
    fit.display = display;
    fit.alpha = 0.1;
    fit.fitExploreSpeed = 1;
    fit.fit2DFirst = 1;
    % }}}

    % Analysis Parameters
    analysis.runAnalysis = runAnalysis;

    % Tracking Parameters
    tracking.runTracking = runTracking;

    params.estimate = estimate;
    params.fit = fit;
    params.tracking = tracking;
    params.analysis = analysis;
    params.LOC = opts.LOC;
    params.CFG = opts.CFG;

end
