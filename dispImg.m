function dispImg(img, order)
        
if isa(img, 'double') || isa(img, 'logical') || isa(img, 'uint8')
    % its a single image
    figure; imagesc(img); axis equal; 
    if size(img,3) == 1
        colormap gray
    end

    pos = get(gcf, 'position');
    set(gcf, 'pos', [-10000 pos(2:end)]);
    set(gcf, 'WindowState', 'maximized');
    set(gca, 'xlim', [1 size(img, 1)], 'ylim', [1 size(img, 2)], 'xtick',[], 'ytick',[] )
    
elseif isa(img, 'cell') && nargin==2
    % there are multiple images for comparison. The display order is
    % specified in the secondary argument 'order' which equals [nrows
    % ncols]
    
    nFig = length(img);
    nR = order(1); nC = order(2);
    figure; 
    
    for jFig = 1 : nFig

        subplot( nR, nC, jFig); imagesc( img{jFig} ); axis equal; 
        if size( img{jFig}, 3) == 1
            colormap gray
        end
        set(gca, 'xlim', [1 size(img{jFig}, 1)], 'ylim', [1 size(img{jFig}, 2)], 'xtick',[], 'ytick',[] )
    end
    pos = get(gcf, 'position');
    set(gcf, 'pos', [-10000 pos(2:end)]);
    set(gcf, 'WindowState', 'maximized');
        
else
    error('Oops, something went wrong!')
end

end