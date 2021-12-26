P00 = [1 0 0 0 0 0 0 1]/2;
P01 = [0 1 1 1 1 1 1 0]/6;
P11 = [0 0 0 1 0 1 1 0]/3;
v01 = P01 - P00;
v11 = P11 - P00;
v10 = v11 - v01;
Pxy = @(x,y) P00 + x*v10 + y*v01;

% compute the plot for the nonnegativity condition
fnonneg = @(x,y) minCoeff(Pxy(x,y));
Inonneg = oracleplot.Implicit.empty([-1.1 1.1], [-0.1 1.1], 'xDivisions', 2^16, 'yDivisions', 2^16);
Inonneg.initializePath(fnonneg, 0, 0.5, 0.05);
Inonneg.closePath(fnonneg);

% compute the plot for the feasible region for the cut inflation
f = @(x,y) slackCut(Pxy(x,y));
I = oracleplot.Implicit.empty([-1.1 1.1], [-0.1 1.1], 'xDivisions', 2^16, 'yDivisions', 2^16);
I.initializePath(f, 0, 0.5, 0.05);
I.closePath(f);

% plot everything
clf;
axis([I.xRange I.yRange]);
hold on;
[x, y] = Inonneg.computePath(f);
fill(x,y,'g');
[x, y] = I.computePath(f);
fill(x,y,'yo');
I.plotCachedPoints;
