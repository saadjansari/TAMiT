%% Monopolar
% get feature
feat = obj.data{jChannel}.features{jTime};
feat.ID = COUNTER;
feat.updateFeatureIDs();
feat.updateFeatureMap();

img = im2double(feat.image);

fs = 12;
% make figure
h = figure('visible', 'on');
hold on; ax = gca; axis ij; set(ax, 'View',[-45,33]); h.Color = 'white';
set(ax,'zlim',[1 size(img,3)])
xlabel('X', 'FontWeight','bold'); ylabel('Y', 'FontWeight','bold'); zlabel('Z', 'FontWeight','bold')
feat.displayFeature3D( ax, size(img,3)); set(ax,'FontSize',fs)
%xlim([60,90]); ylim( [60,100]); zlim([5,7]);
tcks = linspace(0,1, size(img,3)); tcks = tcks(1:2:end);
tcks_label = 1:2: size(img,3);
set(gca,'FontSize',20)
pos = get(gca, 'Position'); set(gca,'Position', [pos(1), pos(2), pos(3), pos(3)*size(img,1)/size(img,2)]);
% cb = colorbar('Location','eastoutside', 'Ticks',tcks, 'TickLabels',tcks_label,  ...
%     'FontWeight', 'bold', 'FontSize', fs, 'Box', 'off', 'LineWidth',1);

%% Bipolar
% get feature
feat = obj.data{jChannel}.features{jTime};
feat.ID = COUNTER;
feat.updateFeatureIDs();
feat.updateFeatureMap();

img = im2double(feat.image);

fs = 12;
% make figure
h = figure('visible', 'on'); 
hold on; ax = gca; axis ij; grid minor; set(ax, 'View',[-45,33]); h.Color = 'white';
set(ax,'zlim',[1 size(img,3)])
xlabel('X', 'FontWeight','bold'); ylabel('Y', 'FontWeight','bold'); zlabel('Z', 'FontWeight','bold')
% set( ax, 'xtick', [], 'ytick', []);
feat.displayFeature3D( ax); set(ax,'FontSize',fs)

% Set All font sizes
for ii = 1:length(h)
    set(h(ii),'FontSize',fs)
end

%% Budding yeast
% get feature
feat = obj.data{jChannel}.features{jTime};
feat.ID = COUNTER;
feat.updateFeatureIDs();
feat.updateFeatureMap();

img = im2double(feat.image);

fs = 12;
% make figure
h = figure('visible', 'on'); 
hold on; ax = gca; axis ij; grid minor; set(ax, 'View',[-45,33]); h.Color = 'white';
set(ax,'zlim',[1 size(img,3)])
xlabel('X', 'FontWeight','bold'); ylabel('Y', 'FontWeight','bold'); zlabel('Z', 'FontWeight','bold')
% set( ax, 'xtick', [], 'ytick', []);
feat.displayFeature3D( ax); set(ax,'FontSize',fs)

% Set All font sizes
for ii = 1:length(h)
    set(h(ii),'FontSize',fs)
end





