function params = initConfiguration( opts)
    % Configuration File: sets up the flags and the paths for running the software in different envoronments

    runFit = 1;
    runAnalysis = 1;
    runTracking = 1;
    
    switch opts.CFG
        case 'RELEASE'
            runLocalFit = 0;
            runGlobalFit = 1;
            runFeatureNumberFit = 1;
            display = 0;

        case 'DEBUG'
            runLocalFit = 1;
            runGlobalFit = 1;
            runFeatureNumberFit = 1;
            display = 1;
    end

    % Estimate Parameters
    estimate.display = display;

    % Spindle
    estimate.spindle.spindleMT = 1;
    estimate.spindle.spindlePoles = 1;
    estimate.spindle.astralMT = 0;
    estimate.spindle.display = display;

    % Monopolar Spindle
    estimate.monopolar.display = display;
    estimate.monopolar.pole = 1;
    estimate.monopolar.mt = 1;

    % Kinetochores
    estimate.kcbank.display = display;

    % Cut7
    estimate.cut7dist.display = display;

    % Fit Parameters
    fit.runFit = runFit;
    fit.runLocalOptimization = runLocalFit;
    fit.runGlobalOptimization = runGlobalFit;
    fit.runFeatureNumberOptimization = runFeatureNumberFit;
    fit.useParallel = false;
    fit.state = opts.CFG;
    fit.display = display;
    fit.alpha = 0.05;
    fit.fitExploreSpeed = 0;

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
