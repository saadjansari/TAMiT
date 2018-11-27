function extractCellsFromMovie( moviePath, whichCells)
% extractCellsFromMovie: extracts information about the cells in a movie.

% Detailed Explanation: This is the first function that is executed in this curved microtubule analysis pipeline. It takes on a multiple of tasks, all associated with pre-processing.
% 1. Identify the movie that will be analyzed in the current run. If moviePath is provided, the identification is complete, and if not provided, the user is prompted in an appropriate directory to select the movie for analysis. At the end of this step, the identification is complete.  
% 2. The movie needs to be imported for processing. A sub-directory is checked for the existence of a previously processed file that has been saved. If preprocessed file exists, no other steps are performed. If no preprocessed file exists, this step takes care of preprocessing. The movie (.nd2 format) is loaded through the bfmatlab package and the data is converted from int16 to int8 for memory conservation. The appropriate variables are created and saved in a .mat file 
% 3. The movie has been loaded but needs to be segmented. If pre-processed file exists and contains segmented cells, this step is skipped. If no pre-processed file exists, a segmentation routine is run and all cells are extracted in 5D (3-space + 1-time + 2-channels). Once segmentation is complete, the segmented cells are stored in the .mat file. The memory is cleared.
% 4. If the input argument whichCells = 'all', all cells are selected for analysis. If whichCells = 'prompt', the user is prompted with a figure with numbered cells and asked to specify cells for analysis.


% Check Inputs
% number of input arguments
if nargin ~= 2
error( sprintf('extractCellsFromMovie: function expects 2 input arguments. input arguments provided = %d', nargin) );
end
% first argument is a string
if ~strcmp( class( moviePath), 'char') 
error( ['extractCellsFromMovie: first argument is expected to be a path of class string, input provided is of type ' class( moviePath) ]);
end
% second  argument is a string matching either 'all' or 'prompt', or 'random'
if ~strcmp( whichCells , 'all') && ~strcmp( whichCells , 'prompt') && ~strcmp( whichCells , 'random') 
error( 'extractCellsFromMovie: second input argument represent which cells to use for analysis. Acceptable strings are "all", "prompt" and "random".')
end
% parse input arguments
if strcmp( moviePath, '')
promptUserForMovie = 1;
else
promptUserForMovie = 0;
end 

startPath = '/Users/saadjansari/Documents/Projects/FY Datasets/';

% Step 1
%  Identify movies
if promptUserForMovie
	
if exist( startUIPath, 'dir') == 7
	[ filename, pathname] = uigetfile( [startUIPath,'*.nd2'], 'Select a .nd2 file');
	else
	[ filename, pathname] = uigetfile( '*.nd2', 'Select a .nd2 file');
	end
	if filename == 0; error('File must be selected for import'); end
	[~, onlyName, ext] = fileparts( filename);
        moviePath = [pathname, filename];
else
	[ pathname, onlyName, ext] = fileparts( moviePath);
end

% Step 2
matPath = [ pathname, filesep, onlyName, '.mat'];
% Check if pre-processed movie exists, if not, load the .nd2 movie and store it in a .mat file
if exist( matPath) ~= 2
	disp('Importing the .nd2 file...')
	% pre process this movie
	importND2( moviePath);
else
	disp('Pre-processed file found in system. Importing will not be performed again.')
end

% Step 3
% We check if movie has been segmented already. We'll probe the saved file from step 2 to check if segmentation has been already performed. If not performed, we'll perform segmentation on the movie
varInfo = who('-file', matPath);
mov = matfile( matPath, 'Writable', true);
if ~ismember('segmentationCheck', varInfo) || mov.segmentationCheck==0 
	% segment this movie
	mov = matfile( matPath, 'Writable', true);
	sizeData = size( mov, 'imData');
	imXYT = squeeze( mov.imData( :, :, 3, 1:5:sizeData(4), 1) );	
	imXY = mean( imXYT, 3);
	disp('Generating the segmentation mask for the full-field movie...')
	maskInfo = generateSegmentationMask( imXY);
	mov.MaskLogical = maskInfo.MaskLogical;
	mov.MaskColored = maskInfo.MaskColored;
	mov.NumCells = maskInfo.NumCells;
	cellsToSegment = [1: maskInfo.NumCells];
	disp('Applying the segmentation mask to the full-field movie...')
	segCells = useSegmentationMask(mov, maskInfo.MaskLogical , cellsToSegment);
	for jCell = 1: length(cellsToSegment)
		currCell = cellsToSegment( jCell);
		mov.( ['cell_' num2str(currCell)]) = segCells(jCell).cell3D;
		mov.( ['cellRaw_' num2str(currCell)]) = segCells(jCell).raw;
	end
	mov.imageNumbered = segCells(1).locations
	mov.segmentationCheck = 1;
	disp(['Segmentation complete. Data stored in location: ', matPath])
else
	disp('Segmented Cells found in system. Segmentation will not be performed again.')
end

% Step 4
if ~ismember('cellCentroids', varInfo)
    imMask = mov.MaskLogical;
    cc = bwconncomp( imMask); st = regionprops( cc, 'Centroid');
    mov.cellCentroids = flip( cat( 1, st.Centroid), 2);
end

end
