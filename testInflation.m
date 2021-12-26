P00 = [1 0 0 0 0 0 0 1]/2;
P01 = [0 1 1 1 1 1 1 0]/6;
P11 = [0 0 0 1 0 1 1 0]/3;
v01 = P01 - P00;
v11 = P11 - P00;
v10 = v11 - v01;
Pxy = @(x,y) P00 + x*v10 + y*v01;

fig = figure;
% compute the plot for the nonnegativity condition
f_nonneg = @(x,y) minCoeff(Pxy(x,y));
data_nonneg = oracleplot.Implicit.empty([-1.1 1.1], [-0.1 1.1], 'xDivisions', 2^16, 'yDivisions', 2^16);
M_nonneg = oracleplot.ImplicitMaker(data_nonneg, f_nonneg, 'figure', fig, 'plotEvaluations', true, 'plotIsolines', true);
M_nonneg.isoline(0, 'featureSize', 0.1);

% compute the plot for the feasible region for the cut inflation
f = @(x,y) slackCut(Pxy(x,y));
data = oracleplot.Implicit.empty([-1.1 1.1], [-0.1 1.1], 'xDivisions', 2^16, 'yDivisions', 2^16);
M = oracleplot.ImplicitMaker(data, f, 'figure', fig, 'plotEvaluations', true, 'plotIsolines', true);
M.isoline(0, 'featureSize', 0.2);
%I.initializePath(f, 0, 0.5, 0.05);
%I.closePath(f);

% plot everything
%clf;
%axis([I.xRange I.yRange]);
%hold on;
%[x, y] = Inonneg.computePath(f);
%fill(x,y,'g');
%[x, y] = I.computePath(f);
%fill(x,y,'yo');
%I.plotCachedPoints;
