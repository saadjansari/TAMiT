[x,y] = meshgrid(10:100,10:100);
z = (1/10)*exp(-((50*(30-x)+50*(-30+y)).^2)/500000).* sqrt(pi/2).*( -erf((5000*(30-x)+5000*(30-y))./(50000*sqrt(2))) + erf((5000*(80-x)+5000*(80-y))./(50000*sqrt(2))) );

figure;
surf(x,y,z); colormap spring;
zticks([]); xticks([]); yticks([]);
xlim([10 100]); ylim([10 100]);
zlabel('Intensity')
xlabel('y')
ylabel('x')
set(gca, 'FontSize', 18);
set(gcf,'Color','w')

%% animation

[x,y] = meshgrid(10:100,10:100);
figure;
z1 = exp(-((1.5*(20-x)+1.5*(-20+y)).^2)/(2*1.5^2 * 100)).* sqrt(pi) .* 10.*10 ./ (2* 1.5*10*sqrt(2));
z2 = -erf((1.5*100*(20-x)+1.5*100*(20-y))./(1000*1.5*sqrt(2))) + erf((1.5*100*(1.5+20-x)+1.5*100*(1.5+20-y))./(1000*1.5*sqrt(2))) ;
z = z1.*z2;
surf(x,y,z); colormap jet; set(gca, 'View', [-71  43.2193])
zticks([]); xticks([]); yticks([]);
xlim([10 100]); ylim([10 100]); zlim([0 max(z(:))])
zlabel('Intensity')
xlabel('y')
ylabel('x')
set(gca, 'FontSize', 18);
set(gcf,'Color','w')
axis tight manual 
set(gca,'nextplot','replacechildren'); 

v = VideoWriter('line2dgaussian.avi');
v.FrameRate = 50;
open(v);

for k = 2:0.5:60
    z1 = exp(-((k*(20-x)+k*(-20+y)).^2)/(2*k^2 * 100)).* sqrt(pi) .* 10.*10 ./ (2* k*10*sqrt(2));
    z2 = -erf((k*100*(20-x)+k*100*(20-y))./(1000*k*sqrt(2))) + erf((k*100*(k+20-x)+k*100*(k+20-y))./(1000*k*sqrt(2))) ;
    z = z1.*z2;
    surf(x,y,z); colormap jet; set(gca, 'View', [-71   43.2193]); zticks([]); xticks([]); yticks([]);
    xlim([10 100]); ylim([10 100]); zlim([0 max(z(:))])
    frame = getframe(gcf);
    writeVideo(v,frame);
    disp(k)
end

close(v);

