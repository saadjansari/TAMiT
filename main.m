function status = main()

    clc; close all;
    clearvars 
    addpath( genpath(pwd) )

    CFG = 'Local';
    settingsFun = 'makeSettings';

    % check if user settings exist in the current folder
    if exist( settingsFun ) ~= 2
        error( ['copy' settingsFunc  '.m to the current folder location : ', pwd] ), end

    % run the settings file : creates a setting.mat file in the results folder
    settingsPath = feval( settingsFun, CFG);

    % Run Single Cell
    singleCell( settingsPath);

    % clear user settings init file from the folder
%     delete [ settingsFun '.m']

end
