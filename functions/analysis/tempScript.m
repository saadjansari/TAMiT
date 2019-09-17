% Create data
nSam = 5;
xstart = randi( 10, [1 nSam]);
xend = randi( [20 100], [1 nSam]);
nPts = randi( [50 70], [1 nSam]);

for jS = 1 : nSam
    xdata{jS} = linspace( xstart(jS), xend(jS), nPts(jS));
    ydata{jS} = xdata{jS} + 10*randn([1 nPts(jS)]);
end

AnalysisSingleCell.plotAreaShadedError(ydata, xdata)
