f = @(x,y) x^2+y^2-0.5;
I = oracleplot.Implicit.empty([-1 1], [-1 1], 'xDivisions', 2^16, 'yDivisions', 2^16);
I.initializePath(f, 0, 0, 0.05);
I.closePath(f);

clf;
axis([I.xRange I.yRange]);
hold on;
[x, y] = I.computePath(f);
fill(x,y,'g');
I.plotCachedPoints;
