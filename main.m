function status = main()

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

    % Run Single Cell
    singleCell( paramsPath);

    % clear user settings init file from the folder
%     delete [ settingsFun '.m']

end
