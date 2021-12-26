f = @(x,y) x^2-y^2;
data = oracleplot.Implicit.empty([-pi pi], [-pi pi], 'xDivisions', 2^16, 'yDivisions', 2^16);

fig = figure;
M = oracleplot.ImplicitMaker(data, f, 'figure', fig, 'plotEvaluations', true, 'plotIsolines', true);
M.isolinesUsingGrid(1, 0.4);
M.isolinesUsingGrid(0.5, 0.4);
M.isolinesUsingGrid(1.5, 0.4);

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
