close all;
% 2,4

N = 6;
load( sprintf('imgInt%.0f.mat',N) );
coords = Methods.FindCurves( imageIn, 'Plot', 1, 'FieldOfView', 60); 
