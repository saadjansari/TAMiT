function status = main( varargin)
    % ----------------------------- PREP ---------------------------
    
    opts = parseArgs( varargin{:} );

    clc; close all;
    clearvars -except CFG opts

    % Make sure we have the correct paths
    warning('off', 'MATLAB:rmpath:DirNotFound');
    rmpath( genpath(pwd) );
    warning('on', 'MATLAB:rmpath:DirNotFound');
    addpath( pwd);

    % check if params file exists in the current folder
    if exist( fullfile( pwd, 'initParams.m' ) ) ~= 2
        error( ['copy initParams.m to the current folder location : ', pwd] );
    end

    % run the settings file : creates a params.mat file in the save directory folder
    paramsPath = feval('initParams', opts);

    % Define cleanup tasks
    c2 = onCleanup( @() delete( gcp('nocreate') ) );
    c3 = onCleanup( @() disp('Closing files and cleaning up') );

    % ----------------------------- MAIN ---------------------------
   
    for jCell = 1 : length( paramsPath)

        % Run Single Cell
        singleCell( paramsPath{ jCell} );

    end

    % ---------------------------- CLEANUP -------------------------
    
    % parseArgs {{{
    function opts = parseArgs( varargin)

        % default
        defaultCFG = 'RELEASE';
        defaultLOC = 'Local';

        % valid Arguments
        validCFG = @(x) strcmpi( x, 'RELEASE') || strcmpi( x, 'DEBUG');
        validLOC = @(x) strcmpi( x, 'Local') || strcmpi( x, 'Summit') || strcmpi( x, 'Rumor');

        % Input Parser
        p = inputParser;
        addParameter( p, 'CFG', defaultCFG);
        addParameter( p, 'LOC', defaultLOC);

        parse( p, varargin{:});
        opts = p.Results;

    end
    % }}}
    
end
