f = @(x,y) x^2+y^2-0.5;
data = oracleplot.Implicit.empty([-1 1], [-1 1], 'xDivisions', 2^16, 'yDivisions', 2^16);
maker = oracleplot.ImplicitMaker(data, f, 'figure', figure, 'plotEvaluations', true);
maker.isoline(0);
maker.isoline(0.1);
maker.isoline(-0.1);