function params = initConfiguration( opts)
    % Configuration File: sets up the flags and the paths for running the software in different environments
    
    % Verbose
    % 0: Essential Logs (Process OFF, Warnings OFF, Messages OFF)
    % 1: Minimal Logs (Process ON, Warnings OFF, Messages OFF)
    % 2: Detailed (Process ON, Warnings ON, Messages ON)
    verbose = 1;
    
    % Use parallel toolbox for fitting
    parallel = 0;
    
    % First fit 2D features (this will lead to more optimization steps but
    % can increase accuracy)
    runFit2DFirst = 1;
    
    % Local Fit (Flag to fit each feature individually)
    runLocalFit = 0;
    
    % Feature Number Optimization (Find statistically optimal number of
    % features by iteratively adding/removing)
    runFeatureNumberFit = 1;
    
    % Visual display
    switch opts.Display
        case 0
            display = 0;
        case 1
            display = 1;
    end
    
    % Common parameters
    common.verbose = verbose;
    common.sigma = [1.2,1.2,1.0];               % Floats - Initial gaussian sigma for features
    
    % Estimation {{{
    % Estimate Parameters
    estimate.display = display;
    estimate.use_time_averaged_image = 1;       % Flag - use a 3 frame time-averaged image for estimation

    % FY Spindle
    estimate.spindle.display = display;             
    estimate.spindle.spindlePoles = 1;                      % Flag - find spindle poles
    estimate.spindle.astralMT = 1;                          % Flag - find astral microtubules
    estimate.spindle.spindleExclusionRange = deg2rad(45);   % Float - cone angle around spindle to avoid searching for astral microtubules 
    estimate.spindle.astralMTRepresentation = 'spherical';  % String - geometric representation (spherical or cartesian)
    estimate.spindle.astralMinLength = 4;                   % Float - minimum length of estimated astral microtubules (keep above 4) 
    estimate.spindle.linewidth = 3;                         % +Int - estimated width of spindle
    estimate.spindle.brightestPixelAsSPB = 0;               % Flag - force one SPB to be at the brightest pixel
    estimate.spindle.visuals = 0;                           % Flag - visuals (set = 0, unless doing specific debugging)
    estimate.spindle.visuals_path = '';                     % Flag - path to save visuals
    estimate.spindle.expectedMAL = 5;                       % +Int - expected major axis length of spindle (5-10 for best results)
    estimate.spindle.minRegionArea = 6;                     % +Int - minimum spindle area in pixels (6-10 for best results)
    estimate.spindle.common = common;
    
    % FY Monopolar Spindle
    estimate.monopolar.display = display;
    estimate.monopolar.astralMT = 1;                        % Flag - find astral microtubules
    estimate.monopolar.astralMTRepresentation = 'cartesian';  % String - geometric representation (spherical or cartesian)
    estimate.monopolar.poleMinIntensity = 0.75;             % Float - minimum intensity of pole relative to brightest image pixel(between 0-1)
    estimate.monopolar.poleDeterminationSensitivity = 0.2;  % Float - extended region sensitivity (best results between 0.2 and 0.4)
    estimate.monopolar.astralMinLength = 4;                   % Float - minimum length of estimated astral microtubules (keep above 4) 
    estimate.monopolar.visuals = 0;                         % Flag - visuals (set = 0, unless doing specific debugging)
    estimate.monopolar.visuals_path = '';                   % Flag - path to save visuals
    estimate.monopolar.common = common;
    
    % BY Spindle
    estimate.spindleBud.display = display;             
    estimate.spindleBud.spindlePoles = 1;                      % Flag - find spindle poles
    estimate.spindleBud.astralMT = 1;                          % Flag - find astral microtubules
    estimate.spindleBud.spindleMinIntensity = 0.65;            % Float - minimum intensity of pole relative to brightest image pixel(between 0-1)
    estimate.spindleBud.spindleDeterminationSensitivity = 0.4; % Float - extended region sensitivity (best results between 0.2 and 0.4)
    estimate.spindleBud.spindleExclusionRange = deg2rad(45);   % Float - cone angle around spindle to avoid searching for astral microtubules 
    estimate.spindleBud.astral_poly_orderXY = 2;               % +Int - polynomial order of curve in X and Y.
    estimate.spindleBud.curvature_model = 'fourier4';          % String - 'fourier1-4'
    estimate.spindleBud.astralMinLength = 4;                   % Float - minimum length of estimated astral microtubules (keep above 4) 
    estimate.spindleBud.linewidth = 2;                         % +Int - estimated width of spindle
    estimate.spindleBud.brightestPixelAsSPB = 0;               % Flag - force one SPB to be at the brightest pixel
    estimate.spindleBud.spb_astral_dist = 0.15;                % Float - distance in microns b/w spb and astral mts
    estimate.spindleBud.voxel_size = [0.05,0.05,0.5];          % Floats - voxel size in microns
    estimate.spindleBud.visuals = 0;                           % Flag - visuals (set = 0, unless doing specific debugging)
    estimate.spindleBud.visuals_path = '';                     % Flag - path to save visuals
    estimate.spindleBud.expectedMAL = 5;                       % +Int - expected major axis length of spindle (5-10 for best results)
    estimate.spindleBud.minRegionArea = 6;                     % +Int - minimum spindle area in pixels (6-10 for best results)
    estimate.spindleBud.common = common;
    estimate.spindleBud.common.sigma = [2.5,2.5,1.2];          % Floats - Initial gaussian sigma for features
    
    % Interphase
%     estimate.interphase.display = display;
%     estimate.interphase.PolyOrder = [3 3 1];
%     estimate.interphase.Visibility = 10;
%     estimate.interphase.FieldOfView = 60;
%     estimate.interphase.StepSize = 5;
%     estimate.interphase.MinLength = 10;
%     estimate.interphase.MaxCurves = 5;
%     estimate.interphase.CostAcceptLink = 10;
%     estimate.interphase.common = common;

    % Kinetochores
    estimate.kcbank.display = display;
    estimate.kcbank.common = common;

    % Cut7
    estimate.cut7dist.display = display;
    estimate.cut7dist.common = common;
    
    % Sid4
    estimate.sid4.display = display;
    estimate.sid4.common = common;

    % Fit {{{
    % Fit Parameters
    fit.pad_xy = 10;            % Number of padding pixels in XY for fitting
    fit.pad_z = 1;              % Number of padding pixels in Z for fitting
    fit.alpha = 0.1;            % Significance level for adding/removing features
    fit.fitExploreSpeed = 0;    % Custom parameter scan speed
    fit.scaleParameters = 1;    % Parameter normalization
    fit.runLocalOptimization = runLocalFit;
    fit.runFeatureNumberOptimization = runFeatureNumberFit;
    fit.useParallel = parallel; 
    fit.state = display;
    fit.display = display;
    fit.fit2DFirst = runFit2DFirst;
    fit.verbose = verbose;
    % }}}

    params.estimate = estimate;
    params.fit = fit;
    params.LOC = opts.LOC;
    params.DisplayState = display;
    params.verbose = verbose;

end
