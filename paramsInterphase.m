function params = paramsInterphase()
    
    % define className (refered to as bank)
    bankName = 'InterphaseMicrotubuleBank';
    
    % function names {{{
    % define class functions to be used for processing
    func.initializeBank = [ bankName ]; 
    func.estimateFeatures = [ 'EstimateMicrotubules'];
    func.fitFeatures = [ 'FitMicrotubules'];
    func.trackFeatures = [ 'trackMicrotubules'];
    func.analyzeFeatures = [ 'analyzeMicrotubules' ];
    func.sim = 'GaussianCurveMaker2D';
    params.func = func;
    params.fitting.fitfunc = 'GaussianCurveMaker2D';
    % }}}

    % Estimation {{{
    % estimation-specific parameters
    params.estimation.visibility = 10;
    params.estimation.fieldOfView = 40;
    params.estimation.stepSize = 5;
    params.estimation.minMTlength = 10;
    params.estimation.maxMTsearchIterations = 5; % each iteration involves searching for 2 microtubules starting from the same point but with opposite orientations.
    params.estimation.cost_accept_link = 10; 
    params.estimation.plotflags.something = 0;
    % }}}

    % Fitting {{{
    % fit logistical info 
    fitDimension = 2; % 2D
    fitPolyOrder = 2; % 2=quadratic(3 parameters), 3=cubic(4 parameters), ...
    fixMTOC = 1; % fit the xyz position of the the nucleation point of the microtubule 
    fitGlobal = 1; % run a global fit after a local fit routine
    fitParallel = 1; % attempt to run fitting in parallel

    % feature appearance info
    sigInit = [1.5, 1.5, 1.0];
    sigRange = [0.2, 0.2, 0.2];
    bkgInit = 0; % will be measured and updated when image is processed
    bkgRange = 0.1;
    ampInit = 0; % will be measured and updated when image is processed
    ampRange = 0.2;
    lengthMTmin = 0.5; % min ratio decrease from estimate
    lengthMTmax = 2.0; % max ratio increase from estimate
    
    % fit toolbox info 
    maxFunEvals = 1e4;
    optTol = 1e-7;
    maxIter = 15;
    tolFun = 1e-6;
    finDifStepSize = 1e-5;
    stepTol = 1e-8;
    displayType = 'iter';

    % Parsing the parameters into a structure {{{
    fitStructInit.dim = fitDimension;
    fitStructInit.polyOrder = fitPolyOrder;
    fitStructInit.fixStartingPoint = fixMTOC;
    fitStructInit.background= bkgInit;
    fitStructUb.background= bkgInit + bkgRange;
    fitStructLb.background= bkgInit - bkgRange;
    fitStructInit.amplitude = ampInit;
    fitStructUb.amplitude = ampInit + ampRange;
    fitStructLb.amplitude = ampInit - ampRange;
    fitStructInit.std = sigInit;
    fitStructUb.std = sigInit + sigRange; 
    fitStructLb.std = sigInit - sigRange;
    params.fitting.initialization.structInit = fitStructInit;
    params.fitting.initialization.structUb = fitStructUb;
    params.fitting.initialization.structLb = fitStructLb;
    params.fitting.initialization.lengthMin = lengthMTmin;
    params.fitting.initialization.lengthMax = lengthMTmax;
    params.fitting.maxFunEvals = maxFunEvals;
    params.fitting.optTol = optTol;
    params.fitting.maxIter = maxIter;
    params.fitting.tolFun = tolFun;
    params.fitting.finDifStepSize = finDifStepSize;
    params.fitting.stepTol = stepTol;
    params.fitting.display = display;
    params.fitting.polyOrder = fitPolyOrder;
    params.fitting.fixMTOC = fixMTOC;
    params.fitting.fitDim = fitDimension;
    params.fitting.global = fitGlobal;
    params.fitting.parallel = fitParallel;
    % }}}
    
    % }}}
    
    % Tracking {{{
    % track-specific parameters
    params.tracking.something = 1;
    
    % }}}
    
    % Analysis {{{
    % analysis-specific parameters
    params.analysis.something = 1;
    
    % }}}

end
