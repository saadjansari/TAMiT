function status = main( CFG)
    % ----------------------------- PREP ---------------------------
    
    clc; close all;
    clearvars -except CFG

    % Make sure we have the correct paths
    warning('off', 'MATLAB:rmpath:DirNotFound');
    rmpath( genpath(pwd) );
    warning('on', 'MATLAB:rmpath:DirNotFound');
    addpath( pwd);

    if nargin == 0
        CFG = 'Local';
    end

    % check if params file exists in the current folder
    initParams = 'initParams';
    if exist( fullfile( pwd, [initParams '.m'] ) ) ~= 2
        error( ['copy ' initParams '.m to the current folder location : ', pwd] );
    end

    % run the settings file : creates a params.mat file in the save directory folder
    paramsPath = feval( initParams, CFG);

    % Define cleanup tasks
    c2 = onCleanup( @() delete( gcp('nocreate') ) );
    c3 = onCleanup( @() disp('Closing files and cleaning up') );

    % ----------------------------- MAIN ---------------------------
   
    for jCell = 1 : length( paramsPath)

        % Run Single Cell
        singleCell( paramsPath{ jCell} );

    end

    % ---------------------------- CLEANUP -------------------------

end
