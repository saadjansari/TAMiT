% This is a testbud for budding yeast cells

%% BY Cell-specification
cellpath = '/Users/saadjansari/Documents/Projects/ImageAnalysis/BY Datasets';
cells = {'1.10.18 37C 1.tif',...
         '1.10.18 37C 2.tif',...
         '1.10.18 37C 4 curved MTs.tif',...
         '1.24.18 25C 3.tif',...
         '1.24.18 37C 6.tif',...
         '2.7.18 37C 2,3 curved MT.tif'};

n_frames = 2;
frames = get_BY_frames_rand( cellpath, cells, n_frames);

%% Nuclear Mask
%[f,ha] = disp_frames( n_frames, length(cells), frames);

masks = find_nuclear_mask( frames);
% [f,ha] = disp_frames( n_frames, length(cells), frames, masks);

%% Estimation
addpath('/Users/saadjansari/Documents/Projects/ImageAnalysis/SingleCell/classes')
addpath(genpath( '/Users/saadjansari/Documents/Projects/ImageAnalysis/SingleCell/functions/external/u-track'))

spindles = find_spindle( frames);
asters = find_astrals( frames, spindles);

%% Add mask-padding
% pad_xy = 5;
% pad_z = 1;
% Add z-padding
% [framed, asters_padded] = padding_on( frames, asters, pad_xy, pad_z);
%[f,ha] = disp_frames( n_frames, length(cells), framed);

%% Test for errors upon exceeding mask

% simImage = test_mask_error( framed{4}, asters_padded{4});
