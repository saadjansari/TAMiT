% Load and Analyze the 'simFitdata.mat file
cd ../..;
addpath(genpath('functions'))

dxyz_new = {};
FP_new = {};
TN_new = {};
lenFP_new = {};
lenTN_new = {};
lenTP_new = {};
foldpath = 'paper/validation/';
fold_names = {'sim_data'};
for ifil = 1 : length(fold_names)
    load( [foldpath, fold_names{ifil}, filesep, 'simFitdata_Monopolar.mat'])
    if isempty( dxyz_new)
        dxyz_new = dxyz;
    else
        dxyz_new = horzcat(dxyz_new , dxyz);
    end
    
    if isempty( FP_new)
        FP_new = FP;
    else
        FP_new = horzcat(FP_new, FP);
    end
    
    if isempty( TN_new)
        TN_new = TN;
    else
        TN_new = horzcat(TN_new, TN);
    end
    
    if isempty( lenFP_new)
        lenFP_new = lenFP;
    else
        lenFP_new = horzcat(lenFP_new, lenFP);
    end
    
    if isempty( lenTN_new)
        lenTN_new = lenTN;
    else
        lenTN_new = horzcat(lenTN_new, lenTN);
    end
    
    if isempty( lenTP_new)
        lenTP_new = lenTP;
    else
        lenTP_new = horzcat(lenTP_new, lenTP);
    end
    
end
dxyz = dxyz_new;
FP = FP_new;
TN = TN_new;
lenFP = lenFP_new;
lenTN = lenTN_new;
lenTP = lenTP_new;

voxel_size = [0.1,0.1,0.5];

% colors
col.rust = [255,112,3]/256;
col.lightblue = [23,141,252]/256;
col.orange = [245,161,29]/256;
col.navy = [4,15,79]/256;
col.yellow = [246,181,27]/256;
col.blue = [18,62,89]/256;
col.cyan = [1,219,218]/256;
col.maroon = [95,0,32]/256;


numNoise = length(snrVals);
data = {};
for nidx = 1 : numNoise
    
    % initialize empty arrays
    err_xyz = [];
    lens_TP = [];
    lens_FP = [];
    lens_FN = [];
    
    for jT = 1 : nTrials
        if ~isempty( dxyz{ nidx, jT} )
            err_xyz = [err_xyz; dxyz{ nidx, jT}.*voxel_size];
        end
        if ~isempty( lenTP{nidx,jT})
            lens_TP = [ lens_TP; sqrt( sum( (lenTP{nidx,jT}.*voxel_size).^2,2))];
        end
        if ~isempty( lenFP{nidx,jT})
            lens_FP = [ lens_FP; sqrt( sum( (lenFP{nidx,jT}.*voxel_size).^2,2))];
        end
        if ~isempty( lenTN{nidx,jT})
            lens_FN = [ lens_FN; sqrt( sum( (lenTN{nidx,jT}.*voxel_size).^2,2))];
        end
    end
    data{nidx}.err_xyz = err_xyz;
    data{nidx}.lens_TP = lens_TP;
    data{nidx}.lens_FP = lens_FP;
    data{nidx}.lens_FN = lens_FN;
    data{nidx}.num_TP = size(lens_TP,1);
    data{nidx}.num_FP = size(lens_FP,1);
    data{nidx}.num_FN = size(lens_FN,1);
end
clearvars -except data snrVals numNoise col voxel_size FP dxyz

%% Remove SNR < 0.6

idxRm = find( snrVals < 0.6);
if ~isempty(idxRm)
    fprintf('Removing SNR values below 0.6...\n')
    
    snrVals(idxRm) = [];
    numNoise = length(snrVals);
    data(idxRm) = [];
    FP(idxRm,:) = [];
    dxyz(idxRm,:) = [];
end

%% Plot: Histograms (XY + Z overlapping)

% bins = linspace(0,0.6,23);
% faceAlpha = 0.7;
% for nidx = 1 : numNoise
% 
%     f = figure('visible', 'on'); f.Color = 'white';
%     % histogram of error_x
%     h1=histogram( (data{nidx}.err_xyz(:,1)+data{nidx}.err_xyz(:,2))/2, bins,  'DisplayName', 'XY');
%     h1.FaceColor = col.maroon; h1.FaceAlpha = faceAlpha; hold on;
%     % histogram of error_z
%     h2=histogram( data{nidx}.err_xyz(:,3), bins,  'DisplayName', 'Z');
%     h2.FaceColor = col.cyan; h2.FaceAlpha = faceAlpha; hold on
%     
%     set(gca,'YScale','log'); ylim([ 0 5000]); xlim([-0.05 0.65]); ylabel('Count'); hold on; 
%     xlabel('Error (microns) ')
%     legend
%     suptitle(['Noise = ' num2str(snrVals(nidx))])
%     set(gca, 'FontSize', 18);
% end

%% Bar plot XY+Z counts
% edges = linspace(0,0.5,13);
% bin_centers = edges(1:end-1)+diff(edges);
% for nidx = 1 : numNoise
%     f = figure('visible', 'on'); f.Color = 'white';
%     
%     [cntXY,~] = histcounts( (data{nidx}.err_xyz(:,1)+data{nidx}.err_xyz(:,2))/2, edges);
%     [cntZ,~] = histcounts( data{nidx}.err_xyz(:,3), edges);
%     % histogram of error_x
%     b = bar(bin_centers, [cntXY;cntZ], 1,'FaceColor', 'flat');
%     set(b(1),'CData',col.yellow, 'DisplayName', ' XY');
%     set(b(2),'CData',col.blue, 'DisplayName', ' Z');
%     
%     set(gca,'YScale','log');  xlim([0.0 0.52]); ylabel('Count'); hold on; 
%     xlabel('Error (\mum) ')
%     legend
%     title(['Noise = ' num2str(snrVals(nidx))])
%     set(gca, 'FontSize', 26);
%     ylim([ 1 5000]);
%     yticks([1,10,100,1000])
% end
%% Plot: Violin plot of XY, and Z errors

noises = [];
err = [];
for nidx = 1 : numNoise
    for idx = 1: size( data{nidx}.err_xyz,1)
        noises = [noises; snrVals( nidx)];
        err = [err; norm( data{nidx}.err_xyz(idx,1:2))];
    end
end
figure;
vs = violinplot(err, noises,...
    'BoxColor',[1,1,1],'EdgeColor',[1,1,1],'ViolinAlpha',0.08, 'ShowData', true, 'Width', 0.4);
xlabel('SNR'); ylabel('Error XY (\mum)')
set(gca, 'FontSize', 24);
set(gcf,'Color','w')
ylim([-0.02,0.62])
yticks([0,0.2, 0.4, 0.6])
xlim([0.5,0.5+ length( unique(noises) )])
%axis tight manual;

noises = [];
err = [];
for nidx = 1 : numNoise
    for idx = 1: size( data{nidx}.err_xyz,1)
        noises = [noises; snrVals( nidx)];
        err = [err; norm( data{nidx}.err_xyz(idx,3))];
    end
end
figure;
vs = violinplot(err, noises,...
    'BoxColor',[1,1,1],'EdgeColor',[1,1,1],'ViolinAlpha',0.08, ...
    'ShowData', true, 'Width', 0.4);
xlabel('SNR'); ylabel('Error Z (\mum)')
set(gca, 'FontSize', 24);
set(gcf,'Color','w')
ylim([-0.02,0.62])
yticks([0,0.2,0.4,0.6])
xlim([0.5,0.5+ length( unique(noises) )])

%% Plot: Violin plot of data

noises = [];
err = [];
for nidx = 2 : numNoise
    for idx = 1: size( data{nidx}.err_xyz,1)
        noises = [noises; snrVals( nidx)];
        err = [err; norm( data{nidx}.err_xyz(idx,:))];
    end
end
figure;
vs = violinplot(err, noises,...
    'BoxColor',[1,1,1],'EdgeColor',[1,1,1],'ViolinAlpha',0.08, 'ShowData', true, 'Width', 0.4);
xlabel('SNR'); ylabel('Mean error (\mum)')
set(gca, 'FontSize', 24);
set(gcf,'Color','w')
ylim([-0.02,0.62])
yticks([0,0.2,0.4,0.6])
xlim([0.5,0.5+ length( unique(noises) )])
%axis tight manual;

%% Histograms: Length of true positives, false positives and false negatives
% bins = linspace(0,5.0,23);
% faceAlpha = 0.7;
% 
% for nidx = 1 : numNoise
% 
%     f = figure('visible', 'on'); f.Color = 'white';
%     % histogram of TP
%     h1=histogram( data{nidx}.lens_TP, bins,  'DisplayName', 'TP');
%     h1.FaceAlpha = faceAlpha; hold on;
%     % histogram of FP
%     h2=histogram( data{nidx}.lens_FP, bins,  'DisplayName', 'FP');
%     h2.FaceAlpha = faceAlpha; hold on
%     % histogram of FN
%     h3=histogram( data{nidx}.lens_FN, bins,  'DisplayName', 'FN');
%     h3.FaceAlpha = faceAlpha; hold on
%     
%     set(gca,'YScale','log'); ylim([ 0 1000]); xlim([-0.05 3.5]); ylabel('Count'); hold on; 
%     xlabel('Length (\mum)')
%     legend
%     suptitle(['Noise = ' num2str(snrVals(nidx))])
%     set(gca, 'FontSize', 18);
% end

%% Plot: Length of true positives, false positives and false negatives
% Mean TP,FP,FN length vs SNR
muTP = zeros(3,numNoise);
muFP = zeros(3,numNoise);
muFN = zeros(3,numNoise);
for nidx = 1 : numNoise
    muTP(:,nidx) = [mean( data{nidx}.lens_TP) , std( data{nidx}.lens_TP)/sqrt(length(data{nidx}.lens_TP) ), length(data{nidx}.lens_TP) ];
    muFP(:,nidx) = [mean( data{nidx}.lens_FP) , std( data{nidx}.lens_FP)/sqrt(length( data{nidx}.lens_FP)), length(data{nidx}.lens_FP) ];
    muFN(:,nidx) = [mean( data{nidx}.lens_FN) , std( data{nidx}.lens_FN)/sqrt(length(data{nidx}.lens_FN )), length(data{nidx}.lens_FN) ];
end

figure;
common_props = {'Linewidth', 0.5, 'LineStyle', 'None', 'Marker', 'o', 'MarkerSize', 16};

eb1 = errorbar(1:numNoise,muTP(1,:), muTP(2,:), 'vertical', common_props{:},...
    'DisplayName', ' Correct', 'Color', col.yellow, 'MarkerFaceColor',col.yellow); hold on
eb2 = errorbar(1:numNoise,muFN(1,:), muFN(2,:), 'vertical', common_props{:},...
    'DisplayName', ' Missed', 'Color', col.blue, 'MarkerFaceColor',col.blue); hold on
% eb3 = errorbar(1:numNoise,muFP(1,:), muFP(2,:), 'vertical', common_props{:},...
%     'DisplayName', ' Spurious', 'Color', col.navy, 'MarkerFaceColor',col.navy); hold on

xticks( 1:length(snrVals))
xticklabels(snrVals(1:end));
xlabel('SNR')
ylabel('Length (\mum)')
% title('Effect of noise on error')
set(gca, 'FontSize', 24);
yticks( [0:1.0:4.0])
set(gcf,'Color','w')
axis tight manual; ylim([-0.2 2.5]); xlim([0.5 numNoise+0.5])
legend('Location', 'best')
hold off

%% Accuracy of Detection (Truth table)
acctable = [];
for nidx = 1 : numNoise
    acctable = [acctable; data{nidx}.num_TP/ (data{nidx}.num_TP + data{nidx}.num_FN), data{nidx}.num_FN/ (data{nidx}.num_TP + data{nidx}.num_FN)];
end
acctable = acctable*100;

figure; bg = bar(acctable, 'stacked');
bg(1).set('DisplayName', 'Correct', 'FaceColor', col.yellow, 'BarWidth',0.5)
bg(2).set('DisplayName', 'Missed', 'FaceColor', col.blue, 'BarWidth',0.5)

xticks( 1:length(snrVals))
xticklabels(snrVals(1:end));
xlabel('SNR')
ylabel('Total Detections (%)')
set(gca, 'FontSize', 24); set(gcf,'Color','w')
axis tight manual; ylim([0 103]); yticks([0:25:100]); xlim([0 length(snrVals)+1])
legend('Location', 'southeast')

%% Plot: Percentage of False Positives vs SNR
fp = zeros(2,numNoise);
for nidx = 1 : numNoise
    fps = [];
    for nt = 1 : size(FP,2)
        % percentage of false positives out of total detections
        if FP{nidx,nt} == 0 && isempty( dxyz{nidx,nt})
            continue
        else
            fps = [fps, 100*FP{nidx,nt}/(FP{nidx,nt}+size( dxyz{nidx,nt},2))];
        end
    end
    fp(:,nidx) = [mean(fps), std(fps)/sqrt(length(fps) )];
end


figure;
common_props = {'Linewidth', 1, 'LineStyle', 'None', 'Marker', 'o', 'MarkerSize', 16};
eb1 = errorbar(1:numNoise,fp(1,:), fp(2,:), 'vertical', common_props{:},...
    'Color', col.blue, 'MarkerFaceColor',col.blue); hold on

xticks( 1:length(snrVals))
xticklabels(snrVals(1:end));
xlabel('SNR')
ylabel('Spurious Detections (%)')
% title('Effect of noise on error')
set(gca, 'FontSize', 24);
yticks( [0:25:100]); 

set(gcf,'Color','w')
axis tight manual; ylim([-5 105]); xlim([0.5 length(snrVals)+0.5])
hold off

%% SNR images
f = figure('visible', 'on'); 
h = tight_subplot(1,4, 0.01, 0.3, 0.1);
snrs = [0.75,1.5, 3 ,5];
[imClean,~,~] = simCell(1000,'Monopolar');
for jj = 1:length(snrs)
    set(f, 'currentaxes', h(jj) );
    
    % Add onise
    noisyy = poissrnd( 5, size(imClean))/5;
    imNoisy = imClean + (1/snrs(jj))*noisyy;
    imagesc( h(jj), max(imNoisy , [], 3) ),
    
    PixelSize = voxel_size(1); 
    Scalebar_length = 1;
    xend = size(imNoisy,2)-4; xstart = xend - Scalebar_length/PixelSize; 
    y0 = size(imNoisy,1)-4;
    % x_location and y_location are wherever you want your scale bar to appear.
    line([xstart, xend],[y0,y0], 'Color','w', 'LineWidth', 4)
                    
    colormap( h(jj), gray); axis equal; xlim( [1 size(imNoisy,1)]); ylim( [1 size(imNoisy,2)]); 
    set( h(jj), 'xtick', [], 'ytick', []);
end
f.Color = 'white';
