function status = main()
    % ----------------------------- PREP ---------------------------
    
    clc; close all;
    clearvars 
    addpath( genpath(pwd) )

    CFG = 'Debug';
    initParams = 'initParams';

    % check if params file exists in the current folder
    if exist( fullfile( pwd, [initParams '.m'] ) ) ~= 2
        error( ['copy ' initParams '.m to the current folder location : ', pwd] );
    end

    % run the settings file : creates a params.mat file in the save directory folder
    paramsPath = feval( initParams, CFG);

    % Define cleanup tasks
    c1 = onCleanup( @() eval('diary off') );
    c2 = onCleanup( @() delete( gcp('nocreate') ) );
    c3 = onCleanup( @() disp('Closing files and cleaning up') );

    % ----------------------------- MAIN ---------------------------
    
    % Run Single Cell
    singleCell( paramsPath);

    % ---------------------------- CLEANUP -------------------------

end
