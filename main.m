function status = main()
    % ----------------------------- PREP ---------------------------
    
    clc; close all;
    clearvars 
    addpath( genpath(pwd) )

    CFG = 'Local';
    initParams = 'initParams';

    % check if user settings exist in the current folder
    if exist( initParams, 'file' ) ~= 2
        error( ['copy ' initParams '.m to the current folder location : ', pwd] );
    end

    % run the settings file : creates a params.mat file in the save directory folder
    paramsPath = feval( initParams, CFG);

    % Define cleanup tasks
    c1 = onCleanup( @() eval('diary off') );
    c2 = onCleanup( @() disp('Closing files and cleaning up') );

    % ----------------------------- MAIN ---------------------------
    
        % Run Single Cell
        singleCell( paramsPath);
        fprintf('MATLAB quit unexpectedly. Closing files and cleaning up...\n')

    % ---------------------------- CLEANUP -------------------------

end
