function MicrotubuleBank = fitInterphaseCurves( MicrotubuleBank, image2D, SavePath)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here







numTubes = length( MicrotubuleBank );
fixStartPoint = 1;

% Prepare microtubules for fitting
for jTube = 1 : numTubes
    
    MicrotubuleBank( jTube).source = image2D;
    MicrotubuleBank( jTube) = prepareForFit( MicrotubuleBank( jTube), 'lsqnonlin' , fixStartPoint);
    
end

% Run Fitting routine
for jTube = 1 : numTubes
    
    MicrotubuleBank( jTube) = fitMTCurves( MicrotubuleBank( jTube) );
    
end

% Discard any dead microtubules
if ~isempty( MicrotubuleBank)
    idxDead = find( [MicrotubuleBank.dead]);
    MicrotubuleBank( idxDead) = [];
    numTubes = length(MicrotubuleBank);
end

% plot all fitted microtubules
plot_overlay = 1;
if plot_overlay
    
    colors2 = distinguishable_colors( 3, {'w', 'k'} );

    % plot vars
    trans = 0.5;
    colInit = [ colors2(1, :), trans];
    colFin = [ colors2(2, :), trans];
    minVal = min( image2D(:) );
    maxVal = max( image2D(:) );
    clim = [minVal, maxVal];
    lw = 5;
    
    figure('Name', 'mt_fit_comparison', 'NumberTitle', 'off');
    subplot(221); imagesc( image2D, clim); colormap gray; axis equal; 
    set(gca, 'FontSize', 20, 'xlim', [1 size(image2D,1)], 'ylim', [1 size(image2D,2)]); title('Interphase Cell')
    
    % display image
    subplot(222); imagesc( image2D, clim); colormap gray; axis equal; hold on;
    for jTube = 1 : numTubes
        
        t = linspace(0,1); % paramateric variable
        xinit = polyval( MicrotubuleBank( jTube).estimatedCoef(1,:), t); % init and final curve coordinates
        yinit = polyval( MicrotubuleBank( jTube).estimatedCoef(2,:), t);
        xfin = polyval( MicrotubuleBank( jTube).fitProps.structFinal.XCoef, t);
        yfin = polyval( MicrotubuleBank( jTube).fitProps.structFinal.YCoef, t);
        
        % plot lines
        plot( xinit, yinit, 'Color', colInit, 'LineWidth', lw); 
        
        % show marker at start point
        plot( xinit(1), yinit(1), '*', 'Color', colors2(3, :), 'MarkerSize', 10, 'LineWidth', lw);
        
    end
    % create legend entry
    clear LH L
    LH(1) = plot(nan, nan, '-', 'LineWidth', lw, 'color', colInit); L{1} = 'Estimated MT';
    LH(2) = plot(nan, nan, '*', 'LineWidth', lw, 'MarkerSize', 10, 'color', colors2(3, :) ); L{2} = 'iMTOC';
    legend(LH, L); hold off; set(gca, 'FontSize', 20, 'xlim', [1 size(image2D,1)], 'ylim', [1 size(image2D,2)]);
    set(gcf, 'pos', get(0, 'ScreenSize') );
    title('MT Estimates')
    
    
    
    % display image
    subplot(224); imagesc( image2D, clim); colormap gray; axis equal; hold on;
    for jTube = 1 : numTubes
        
        t = linspace(0,1); % paramateric variable
        % init and final curve coordinates
        xfin = polyval( MicrotubuleBank( jTube).fitProps.structFinal.XCoef, t);
        yfin = polyval( MicrotubuleBank( jTube).fitProps.structFinal.YCoef, t);
        
        % plot lines
        plot( xfin, yfin, 'Color', colFin, 'LineWidth', lw);
        
        % show marker at start point
        plot( xfin(1), yfin(1), '*', 'Color', colors2(3, :), 'MarkerSize', 10, 'LineWidth', lw);
        
    end
    % create legend entry
    clear LH L
    LH(1) = plot(nan, nan, '-', 'LineWidth', lw, 'color', colFin); L{1} = 'Fitted Line';
    LH(2) = plot(nan, nan, '*', 'LineWidth', lw, 'MarkerSize', 10, 'color', colors2(3, :) ); L{2} = 'iMTOC';
    legend(LH, L); hold off; set(gca, 'FontSize', 20, 'xlim', [1 size(image2D,1)], 'ylim', [1 size(image2D,2)]);
    set(gcf, 'pos', get(0, 'ScreenSize') );
    title('MT Fits')
    
    
    
    
    % display image
    subplot(223); imagesc( image2D, clim); colormap gray; axis equal; hold on;
    for jTube = 1 : numTubes
        
        t = linspace(0,1); % paramateric variable
        xinit = polyval( MicrotubuleBank( jTube).estimatedCoef(1,:), t); % init and final curve coordinates
        yinit = polyval( MicrotubuleBank( jTube).estimatedCoef(2,:), t);
        xfin = polyval( MicrotubuleBank( jTube).fitProps.structFinal.XCoef, t);
        yfin = polyval( MicrotubuleBank( jTube).fitProps.structFinal.YCoef, t);
        
        % plot lines
        plot( xinit, yinit, 'Color', colInit, 'LineWidth', lw); 
        plot( xfin, yfin, 'Color', colFin, 'LineWidth', lw);
        
        % show marker at start point
        plot( xinit(1), yinit(1), '*', 'Color', colors2(3, :), 'MarkerSize', 10, 'LineWidth', lw);
        plot( xfin(1), yfin(1), '*', 'Color', colors2(3, :), 'MarkerSize', 10, 'LineWidth', lw);
        
    end
    % create legend entry
    clear LH L
    LH(1) = plot(nan, nan, '-', 'LineWidth', lw, 'color', colInit); L{1} = 'Estimated Line';
    LH(2) = plot(nan, nan, '-', 'LineWidth', lw, 'color', colFin); L{2} = 'Fitted Line';
    LH(3) = plot(nan, nan, '*', 'LineWidth', lw, 'MarkerSize', 10, 'color', colors2(3, :) ); L{3} = 'iMTOC';
    legend(LH, L); hold off; set(gca, 'FontSize', 20, 'xlim', [1 size(image2D,1)], 'ylim', [1 size(image2D,2)]);
    set(gcf, 'pos', get(0, 'ScreenSize') );
    title('Cubic Fit Overlay')
    
end

SaveOpenFigsToFolder(SavePath, 'png', 1);

end

