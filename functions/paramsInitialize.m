function params = paramsInitialize()
% This is a script that initializes the parameters for various processes including estimation, fitting, tracking and plotting

% Cell Phase Classification
cellPhaseClassification.mitosisSNR = 3;
cellPhaseClassification.mitosisMinPixels = 5;

% General Estimation
estimation.InitImageFrameRange = 3;
estimation.usePreviousFitAsEstimate = 1;

% General Fitting
fitting.doFitting = 0;

% General Tracking
tracking.doTracking = 0;

% Interphase (Curved Microtubules) :
paramsIntMT = 'paramsInterphase';
params.interphase = paramsMicrotubulesInterphase(); 

paramsMitMT = 'paramsMicrotubulesMitosis';
params.mitosis = paramsMicrotubulesMitosis();

paramsMitKC = 'paramsKinetochoresMitosis';
params.kc = paramsKinetochoresMitosis();


params.cellPhaseClassification = cellPhaseClassification;
params.estimation = estimation;
params.fitting = fitting;
params.tracking = tracking;

% toggle all plot flags
doToggle = 0;
flagValue = 0;
if doToggle
    [params.interphase.plotflag, params.mitosis.plotflag, params.kc.plotflag] = togglePlotFlags( params.interphase.plotflag, params.mitosis.plotflag, params.kc.plotflag, flagValue);
end

    % paramsMicrotubulesMitosis {{{
    function params = paramsMicrotubulesMitosis()
        
        % define className (refered to as bank)
        bankName = 'MitosisMicrotubuleBank';
       
        % define class functions to be used for processing
        func.initializeBank = [ bankName ]; 
        func.estimateFeatures = [ 'EstimateMicrotubules'];
        func.fitFeatures = [ 'FitMicrotubules'];
        func.trackFeatures = [ 'trackMicrotubules'];
        func.analyzeFeatures = [ 'analyzeMicrotubules' ];
        func.sim = 'GaussianCurveMaker2D';
        params.func = func;
        
        % estimation-specific parameters
        params.estimation.something = 8;

        % fit-specific parameters
        params.fitting.something = 1;


        % track-specific parameters
        params.tracking.something = 1;


        % analysis-specific parameters
        params.analysis.something = 1;


        % display-specific parameters
        params.plotflag.something = 1;
        params.display.something = 1;
    
    end
    %  }}}

    % paramsKinetochoresMitosis {{{
    function params = paramsKinetochoresMitosis()
        
        % define className (refered to as bank)
        bankName = 'MitosisKinetochoresBank';
       
        % define class functions to be used for processing
        func.initializeBank = [ bankName ]; 
        func.estimateFeatures = [ 'EstimateMicrotubules'];
        func.fitFeatures = [ 'FitMicrotubules'];
        func.trackFeatures = [ 'trackMicrotubules'];
        func.analyzeFeatures = [ 'analyzeMicrotubules' ];
        func.sim = 'GaussianCurveMaker2D';
        params.func = func;
        
        % estimation-specific parameters
        params.estimation.something = 8;

        % fit-specific parameters
        params.fitting.something = 1;


        % track-specific parameters
        params.tracking.something = 1;


        % analysis-specific parameters
        params.analysis.something = 1;


        % display-specific parameters
        params.plotflag.something = 1;
        params.display.something = 1;
    
    end
    %  }}}

    % togglePlotFlags {{{
    function varargout = togglePlotFlags( varargin)
        toggle = varargin{end};
        for jVar = 1: length( varargin)-1
            par = varargin{ jVar};
            varargout{jVar} = par;
            fnames = fieldnames( par);
            for jField = 1 : length(fnames)
                varargout{jVar}.(fnames{jField}) = toggle;
            end
        end
    end
    % }}}

%  ------ The End ------
end
