function showimage(fign, img, text)
% [N, D, M] = GetDimensions(x)
%
% Display an image in a given figure, with a given title.
%
% fign:   The figure number.
% img:    The image.
% img:    The title text.
%
% (c) Cedric Vonesch, Biomedical Imaging Group, EPFL, 2009.07.01

figure(fign);
clf;
image(img);
colormap(gray(256));
axis image;
axis off;
title(text);