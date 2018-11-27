function [ imMicrotubules2D, metaData, SavePath] = quickLoad( method)
% quickLoad: Used to quickly load data from a stored .mat file or run
% loadND2file.m.
% if method = 1, load a stored .mat file in a specific directory from
% datasets. This file would be a result of having run loadND2file.m earlier
% on and storing the useful result.
% if method = 2, run loadND2file.m


% TO-DO: pre-process all the movies and save them in the appropriate format. When running the detection, load from the pre-processed datastore if the file is present in the store, otherwise run the processing function, add the movie information to the data storage and continue.

% Ideas: 
% For better efficiency save the segmented data as well, and only load it when it is needed. But first i need to be happy with the segmentation procedure. happy is a strong word. lets say i need to be reasonably satisfied. Am i reasonably satisfied? I think its important to store it because segmented 3D movies can practically take up as much space as a 3D 512x512 movie. I need some way of storing that data. Since i process things cell by cell (and frame by frame), it makes sense to load up movies cell by cell, rather than loading up the entire list of cells. This will be crucial when i extrapolate to using 3D information.

if method == 1
    
    % set directory to directory containing stored .mat file
    cDir = '/Users/saadjansari/Documents/Projects/FY Datasets/quickLoad/';
    load([cDir, 'temp.mat'])
    
elseif method == 2
    
    % load the .nd2 file
    fileInfo = loadND2file();
    
    % process useful channel and save useful information
    imMicrotubules = fileInfo.img( :,:,:,:,1); % 1st channel is microtubules
    imMicrotubules2D = squeeze( max( imMicrotubules, [], 3) ); % take a max intensity projection in the z-dimension
    metaData = fileInfo.metaData;
    
    % Now we have a 3-dimensional array with (x,y,t)
    % save it for temporary loading later
    cDir = '/Users/saadjansari/Documents/Projects/FY Datasets/quickLoad/';
    if exist(cDir, 'dir') ~= 7
        mkdir( cDir);
    end
    save([cDir, 'temp'], 'imMicrotubules2D', 'metaData'); 
    
end

movieName = metaData.fileName;

folderTimeStamp = datestr( now, 'yymmdd_HHMM');
SavePath = [ pwd, filesep, 'Results', filesep, folderTimeStamp];

SavePath = [ SavePath, filesep, movieName];

end
