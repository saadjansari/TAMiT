function [h,ha] = disp_frames( nRows, nCols, frames, masks)
    
    if nRows * nCols < length( frames)
        error('nRows*nCols must be greater than the number of frames')
    end
    addpath('/Users/saadjansari/Documents/Projects/ImageAnalysis/SingleCell/functions/external/tight_subplot')
    h = figure;
    ha = tight_subplot( nRows, nCols, 0.01, 0.01, 0.01);
    
    % Apply mask if available
    masked_frames2D = frames;
    if nargin == 4
        for jframe = 1 : length(masked_frames2D)
            masked_frames2D{jframe} = mat2gray(frames{jframe}) .* masks{jframe};
        end
    end
    
    % Display frames
    for jframe = 1 : length(masked_frames2D)
        
        % 2D image
        im2 = mean( masked_frames2D{jframe}, 3);
        axes( ha( jframe) );
        imagesc( im2 );
        colormap gray
        axis equal; 
        xlim( [1 size(im2, 2) ]); 
        ylim( [1 size(im2, 1) ]); 
        set( ha( jframe), 'xtick', [], 'ytick', []);
    end
end