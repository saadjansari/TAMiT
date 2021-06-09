function make_fig6_curvature( img, feats, pth)

    ttt = [22,23,24,25,26,27];
    timeStep = 5;
    pix_size = [0.05,0.05];
    custom_xy = 1;
    
    f = figure('visible', 'on');
    h = tight_subplot(2,length(ttt), [-0.32,0.002]);
    fs = 16;
    set(f,'Position', get(groot,'Screensize'))

    for jTime = 1:length(ttt)
        
        ct = ttt(jTime);
        
        % Image to display
        im2 = img(:,:,ct);
        imgauss = imgaussfilt(im2,1);
        med = median(imgauss(im2(:)~=0));
        J = imadjust(imgauss, [med, 4*med],[]);
        
        % Custom XY range
        xrange = 1: size(img(:,:,1),1); yrange = 1:size(img(:,:,1),2);
        if custom_xy
            yrange = 15:80; xrange = 20:85;
        end
        
        % Axes 1 : Original image
        axn = jTime;
        set(f, 'currentaxes', h(axn) );
        
        imagesc( h(axn), J )
        colormap( h(axn), gray); axis equal; xlim( [xrange(1) xrange(end) ]); ylim( [yrange(1) yrange(end) ]); 
        set( h(axn), 'xtick', [], 'ytick', []);
        title({['Frame = ' num2str( ct)],['Time = ',num2str( (ttt(jTime)-ttt(1))* timeStep),' sec']})
        hold on;

        % Scale bar
        if axn == 1
            PixelSize = pix_size(1); Scalebar_length = 1;
            xend = xrange(end)-4; xstart = xend - Scalebar_length/PixelSize; y0 = yrange(end)-4;
            % x_location and y_location are wherever you want your scale bar to appear.
            line([xstart, xend],[y0,y0], 'Color','w', 'LineWidth', 4)
        end
        
        % Axes 2: overlayed features
        axn = jTime+length(ttt);
        set(f, 'currentaxes', h(axn) );
        
        % Display image
        imagesc( h(axn), J )
        colormap( h(axn), gray); axis equal; xlim( [xrange(1) xrange(end) ]); ylim( [yrange(1) yrange(end) ]); 
        set( h(axn), 'xtick', [], 'ytick', []);
        hold on;

        % Display features
        
        % Display the spindles
        spindles = feats{1};
        spindles{ct}.displayFeature( h(axn));
        % change spindle color/width
        objs = findall( h(axn)); sp_line = objs(2);
        set(sp_line, 'LineWidth',6, 'Color',[0 19 222]/255);
        % Add SPBs
        spb1 = spindles{ct}.startPosition(1:2);
        spb2 = spindles{ct}.endPosition(1:2);
        plot( spb1(1), spb1(2), '.','Color',[0 19 222]/255, 'MarkerSize', 30);
        plot( spb2(1), spb2(2), '.','Color',[0 19 222]/255, 'MarkerSize', 30);
        
        % Display features at each SPB
        bud1 = feats{2};
        bud2 = feats{3};
        bud1.drawMatchedFeature_curvature( h(axn), ct);                            
        bud2.drawMatchedFeature_curvature( h(axn), ct);
        
        % Colorbar
        if ct == ttt(end)
            
            h_temp2 = axes('Position',get(h(axn),'Position'), 'Color', 'none', 'xtick', [], 'ytick', [], 'XColor', 'none', 'YColor', 'none');
            axis equal; axis ij; xlim( [xrange(1) xrange(end) ]); ylim( [yrange(1) yrange(end) ]);
            set(f, 'currentaxes', h_temp2 );hold on
            colormap( h_temp2, brighten(jet,0.4) );

            tcks = linspace(0,1,4);
            tcks_label = linspace(0,1.5,4);
            cb = colorbar('Location','eastoutside', 'Ticks',tcks, ...
                'TickLabels', tcks_label, ...
                'FontSize', fs, 'Box', 'off', 'LineWidth',0.5);
            cb.Label.String = 'Curvature (\mum^{-1})'; 
%             cb.Label.Rotation = 0;
            set(h_temp2, 'Position',get(h(axn),'Position'))

            % Change cb height to be within 90% of axes
            pos0 = get(cb, 'Position'); pos1=pos0;
            cen = (pos0(2) + pos0(2)+pos0(4))/2;
            bottom = cen - 1.01*(cen-pos0(2)); top = cen + 1.01*(cen-pos0(2));
            pos1([2,4]) = [ bottom, top-bottom];
            pos1(3) = 0.01;
            set(cb, 'Position', pos1)
        end
            
%         title(['Time = ', num2str( ct)]);
        drawnow;
        pause(0.1);

    end
    
    % Set All font sizes
    for ii = 1:length(h)
        set(h(ii),'FontSize',fs)
    end
    f.Color = 'white';

end