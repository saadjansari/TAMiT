addpath(genpath('classes'))
invert=1;

%% Spot
img = zeros(100);
myspot = Spot([51,51],10,[15,15],2,{},{});
imSim = myspot.simulateFeature( size(img));
imSim = imcomplement(imSim);
figure; 
if invert
    imagesc( imcomplement(imSim));
else
    imagesc(imSim);
end
imagesc(imSim); colormap gray; axis equal;
xticks([]); yticks([]); xlim( [1, size(img,1)]); ylim([1, size(img,2)]);


%% Line
img = zeros(100);
myline = Line([89,31],[11,69],10,[10,10],2,{},{});
imSim = myline.simulateFeature( size(img));
figure; 
if invert
    imagesc( imcomplement(imSim));
else
    imagesc(imSim);
end
colormap gray; axis equal;
xticks([]); yticks([]); xlim( [1, size(img,1)]); ylim([1, size(img,2)]);


%% Curve
img = zeros(100);
mycurve = CurvedMT([85,20],0.9*pi,-0.01,93, 50, [5,5],2,{},{});
mycurve.params.idx = [1:numel(img)]'; 
[mycurve.params.x,mycurve.params.y] = ind2sub(size(img), mycurve.params.idx);
imSim = mycurve.simulateFeature( size(img));
figure; 
if invert
    imagesc( imcomplement(imSim));
else
    imagesc(imSim);
end
colormap gray; axis equal;
xticks([]); yticks([]); xlim( [1, size(img,1)]); ylim([1, size(img,2)]);


%% Monopolar spindle
img = zeros(100);
myspot = Spot([51,51],10,[6,6],2,{},{});
mylines{1} = Line([51,51],[91,40],10,[5,5],2,{},{});
mylines{2} = Line([51,51],[15,81],10,[5,5],2,{},{});
mylines{3} = Line([51,51],[21,11],10,[5,5],2,{},{});

myaster = AsterMT( 2, myspot, mylines{:});
imSim = myaster.simulateFeature( size(img));
% imSim = imcomplement(imSim);

figure; 
if invert
    imagesc( imcomplement(imSim));
else
    imagesc(imSim);
end
colormap gray; axis equal;
xticks([]); yticks([]); xlim( [1, size(img,1)]); ylim([1, size(img,2)]);

%% Bipolar spindle
img = zeros(100);
idx = [1:numel(img)]'; [x,y] = ind2sub(size(img), idx);
params = struct('idx',idx, 'x', x,'y',y);

% spindle
spb1 = Spot([40,45],20,[3,3],2,{},{});
spb2 = Spot([65,57],20,[3,3],2,{},{});
myspindle = Line(spb1.position,spb2.position,25,[3,3],2,{},{});

% aster 1
mycurves1 = [];
mycurves1{1} = CurvedMT(flip(spb1.position),1.35*pi,0.06,35, 10, [2,2],2,{},{});
mycurves1{1}.params = params;
mycurves1{2} = CurvedMT(flip(spb1.position),1.3*pi,-0.025,35, 10, [2,2],2,{},{});
mycurves1{2}.params = params;

myaster1 = Aster( 2, spb1, mycurves1{:});


% aster 2
mycurves2 = [];
mycurves2{1} = CurvedMT(flip(spb2.position),0.3*pi,0.04,30, 10, [2,2],2,{},{});
mycurves2{1}.params = params;

myaster2 = Aster( 2, spb2, mycurves2{:});

% bipolar spindle
my_bipolarspindle = SpindleNew( 2, img, {myspindle,myaster1,myaster2}, {});


imSim = my_bipolarspindle.simulateFeature( size(img));
figure; 
if invert
    imagesc( imcomplement(imSim));
else
    imagesc(imSim);
end
colormap gray; axis equal;
xticks([]); yticks([]); xlim( [1, size(img,1)]); ylim([1, size(img,2)]);
