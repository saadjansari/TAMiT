% Simulate image, add noise, then fit, and make pretty pictures

%% Get good image for paper
snrs = [0.75,1.0,1.25,1.5,3.0];
nsnr = length(snrs);

% sim image
[imSim, mobj, params] = simMonopolar();

figure; imagesc( max( imSim,[],3)); colormap gray; axis equal;
title(['SNR = ', num2str(1)]);


%% Add snr to all images and view
figure;
ha = tight_subplot( 1,nsnr,[.01 .03],[.1 .01],[.01 .01]);

for ii = 1: nsnr
    
    imNoisy = max( mat2gray( add_poisson_noise(imSim, snrs(ii)) ),[],3);
    axes(ha(ii) );
    imagesc( imNoisy);
    colormap gray; axis equal;
    xlim([1,size(imSim,1)]); ylim([1,size(imSim,1)]);
    xticks([]); yticks([]);
    title(['SNR = ', num2str( snrs(ii))]);
end
set(ha(1:end),'XTickLabel',''); set(ha(1:end),'YTickLabel','')

%% Fit all
finalFits = cell(nsnr);
imNoisys = cell(nsnr);

for ii = 1: nsnr
    imNoisys{ii} = add_poisson_noise(imSim, snrs(ii) );
    finalFits{ii} = fitImageWithNoise(imSim, mobj, params, snrs(ii) );  
end

%% Display


f = figure; f.Color = 'white'; set(f,'Position', get(groot,'Screensize'))
ha = tight_subplot( 2,nsnr,[-.28 .005],[.0 .0],[.05 .05]);

for ii = 1 : nsnr
    jj = nsnr+ii;
    
    im2 = mat2gray(max( imNoisys{ii},[],3));
    imgauss = mat2gray( imgaussfilt(im2,1) );
    med = median(imgauss(im2(:)~=0));
    J = imadjust(imgauss, [med, 3*med],[]);
    % Original image
    axes(ha(ii) );
    imagesc( im2 );
    colormap gray; axis equal;
    xlim([1,size(imSim,1)]); ylim([1,size(imSim,1)]);
    xticks([]); yticks([]);
    title(['SNR = ', num2str( snrs(ii))]);
    
    % Scale bar
    if ii == 1
        PixelSize = 0.1067; Scalebar_length = 1;
        xend = size(im2,1)-5; xstart = xend - Scalebar_length/PixelSize; y0 = size(im2,2)-5;
        % x_location and y_location are wherever you want your scale bar to appear.
        line([xstart, xend],[y0,y0], 'Color','w', 'LineWidth', 4)
    end
    
    % Features 3D
    axes(ha(jj) );
    imagesc( im2); hold on;
    colormap gray; axis equal;
    xlim([1,size(imSim,1)]); ylim([1,size(imSim,1)]);
    xticks([]); yticks([]);
    
    % Axes 2 for 3D features
    h_temp = axes('Position',get(ha(jj),'Position')); % make new axes
    set(h_temp, 'Color', 'none'); % i thought this may make the new axes background transparent, but it doesn't work 

    axis equal; axis ij;
    set( h_temp, 'xtick', [], 'ytick', []);
    xlim([1,size(imSim,1)]); ylim([1,size(imSim,1)]);
    set(f, 'currentaxes', h_temp );hold on
    colormap( h_temp, hsv);
    finalFits{ii}.displayFeature( h_temp,1);
end
% Colorbar
% Axes 3 for colorbar
h_temp2 = axes('Position',get(ha(end),'Position'), 'Color', 'none', 'xtick', [], 'ytick', [], 'XColor', 'none', 'YColor', 'none');
axis equal; axis ij; xlim([1,size(imSim,1)]); ylim([1,size(imSim,1)]);
set(f, 'currentaxes', h_temp2 );hold on
colormap( h_temp2, hsv);
fs=12;
tcks = linspace(0,1, size(imSim,3)); tcks = tcks(1:2:end);
tcks_label = ([1:2: size(imSim,3)]-1)*0.5;
cb = colorbar('Location','eastoutside', 'Ticks',tcks, 'TickLabels',tcks_label,  ...
    'FontSize', fs, 'Box', 'off', 'LineWidth',0.5);
% cb.Label.String = 'Z(\mum)'; cb.Label.Rotation = 0;
set(h_temp2, 'Position',get(ha(end),'Position'))

% Change cb height to be within 90% of axes
pos0 = get(cb, 'Position'); pos1=pos0;
cen = (pos0(2) + pos0(2)+pos0(4))/2;
bottom = cen - 0.95*(cen-pos0(2)); top = cen + 0.95*(cen-pos0(2));
pos1([2,4]) = [ bottom, top-bottom];
pos1(3) = 0.015;
set(cb, 'Position', pos1)


%% Display (inverted)

f = figure; f.Color = 'white'; set(f,'Position', get(groot,'Screensize'))
ha = tight_subplot( 2,nsnr,[-.28 .005],[.0 .0],[.05 .05]);

for ii = 1 : nsnr
    jj = nsnr+ii;
    
    im2 = mat2gray(max( imNoisys{ii},[],3));
    imgauss = mat2gray( imgaussfilt(im2,1) );
    med = median(imgauss(im2(:)~=0));
    J = imadjust(imgauss, [med, 3*med],[]);
    % Original image
    axes(ha(ii) );
    imagesc( imcomplement(im2) );
    colormap gray; axis equal;
    xlim([1,size(imSim,1)]); ylim([1,size(imSim,1)]);
    xticks([]); yticks([]);
    title(['SNR = ', num2str( snrs(ii))]);
    
    % Scale bar
    if ii == 1
        PixelSize = 0.1067; Scalebar_length = 1;
        xend = size(im2,1)-5; xstart = xend - Scalebar_length/PixelSize; y0 = size(im2,2)-5;
        % x_location and y_location are wherever you want your scale bar to appear.
        line([xstart, xend],[y0,y0], 'Color','k', 'LineWidth', 4)
    end
    
    % Features 3D
    axes(ha(jj) );
    imagesc( imcomplement(im2) ); hold on;
    colormap gray; axis equal;
    xlim([1,size(imSim,1)]); ylim([1,size(imSim,1)]);
    xticks([]); yticks([]);
    
    % Axes 2 for 3D features
    h_temp = axes('Position',get(ha(jj),'Position')); % make new axes
    set(h_temp, 'Color', 'none'); % i thought this may make the new axes background transparent, but it doesn't work 

    axis equal; axis ij;
    set( h_temp, 'xtick', [], 'ytick', []);
    xlim([1,size(imSim,1)]); ylim([1,size(imSim,1)]);
    set(f, 'currentaxes', h_temp );hold on
    colormap( h_temp, hsv);
    finalFits{ii}.displayFeature( h_temp,1);
end
% Colorbar
% Axes 3 for colorbar
h_temp2 = axes('Position',get(ha(end),'Position'), 'Color', 'none', 'xtick', [], 'ytick', [], 'XColor', 'none', 'YColor', 'none');
axis equal; axis ij; xlim([1,size(imSim,1)]); ylim([1,size(imSim,1)]);
set(f, 'currentaxes', h_temp2 );hold on
colormap( h_temp2, hsv);
fs=12;
tcks = linspace(0,1, size(imSim,3)); tcks = tcks(1:2:end);
tcks_label = ([1:2: size(imSim,3)]-1)*0.5;
cb = colorbar('Location','eastoutside', 'Ticks',tcks, 'TickLabels',tcks_label,  ...
    'FontSize', fs, 'Box', 'off', 'LineWidth',0.5);
% cb.Label.String = 'Z(\mum)'; cb.Label.Rotation = 0;
set(h_temp2, 'Position',get(ha(end),'Position'))

% Change cb height to be within 90% of axes
pos0 = get(cb, 'Position'); pos1=pos0;
cen = (pos0(2) + pos0(2)+pos0(4))/2;
bottom = cen - 0.95*(cen-pos0(2)); top = cen + 0.95*(cen-pos0(2));
pos1([2,4]) = [ bottom, top-bottom];
pos1(3) = 0.015;
set(cb, 'Position', pos1)

