classdef ImplicitMaker < handle
% Class containing the plotting logic for `.Implicit`
%
% Works using a variant of the marching squares algorithm.

    properties (SetAccess = protected)
        data % (`.Implicit`): Implicit plot data
        f % (function_handle): Function to plot
        evaluationFunction % (function_handle): Function to use for evaluation, may include plotting
        figure % (handle or ``[]``): Figure handle to display progress, optional
        figureIsolines % (cell(1,\*) of handles): Isoline line handles
        plotEvaluations % (logical): Whether to plot new function evaluations
        plotIsolines % (logical): Whether to plot new isolines
    end

    methods

        function self = ImplicitMaker(data, f, varargin)
        % Creates a plotter object
        %
        % Args:
        %   data (`.Implicit`): Implicit plot data
        %   f (function_handle): Function to plot
        %
        % Keyword Args:
        %   figure (handle): Figure handle to display progress, optional
        %   plotEvaluations (logical): Whether to plot function evaluations, default: false
        %   plotIsolines (logical): Whether to plot isolines, default: true
            args = struct('figure', {[]}, 'plotEvaluations', {false}, 'plotIsolines', {true});
            args = oracleplot.populateStruct(args, varargin);
            self.data = data;
            self.f = f;
            self.figure = args.figure;
            if ~isempty(self.figure)
                figure(self.figure);
                hold on;
            end
            self.plotEvaluations = args.plotEvaluations;
            self.plotIsolines = args.plotIsolines;
            if self.plotEvaluations
                self.evaluationFunction = @(x, y) self.eval(x, y);
            else
                self.evaluationFunction = f;
            end
        end

        function v = eval(self, x, y)
            v = self.f(x, y);
            if ~isempty(self.figure) && self.plotEvaluations
                figure(self.figure);
                plot(x, y, 'gx');
            end
        end
    end

    methods % Isoline creation

% $$$         function indices = isolinesUsingGrid(self, altitude, cellSize)
% $$$         % Finds isolines at a given altitude using a discretization scheme
% $$$         %
% $$$         % Args:
% $$$         %   altitude (double): Altitude at which to compute the isolines
% $$$         %   cellSize (double): Size of the grid cells; should be roughly half the size of the smallest component to find
% $$$         %
% $$$         % Returns:
% $$$         %   integer(1,\*): Indices of the created isolines
% $$$             f = self.evaluationFunction;
% $$$             delta = self.powerOfTwoDelta(cellSize);
% $$$             dX = self.data.xDivisions/delta;
% $$$             dY = self.data.yDivisions/delta;
% $$$             values = zeros(dX+1,dY+1);
% $$$             for i = 0:dX
% $$$                 for j = 0:dY
% $$$                     values(i+1,j+1) = f(i*delta, j*delta);
% $$$                 end
% $$$             end
% $$$
% $$$         end

        function j = isoline(self, altitude, varargin)
        % Creates an isoline at the given altitude
        %
        % Args:
        %   altitude (double): Altitude of the isoline to plot
        %
        % Keyword Args:
        %   below (double(1,2) or ``[]``): Coordinates of a point below the given altitude, optional
        %   above (double(1,2) or ``[]``): Coordinates of a point above the given altitude, optional
        %   featureSize (double): Estimate of the size of the features of the isoline to plot
            args = struct('below', {[]}, 'above', {[]}, 'featureSize', {[]});
            args = oracleplot.populateStruct(args, varargin);
            if ~isempty(args.featureSize)
                deltax = args.featureSize/abs(self.data.xRange(2) - self.data.xRange(1))*self.data.xDivisions;
                deltay = args.featureSize/abs(self.data.yRange(2) - self.data.yRange(1))*self.data.yDivisions;
                delta = 2^floor(log2(min(deltax, deltay)));
            else
                delta = max(min(self.data.xDivisions, self.data.yDivisions)/128, 1);
            end
            if ~isempty(args.below)
                below = args.below;
                [ix, iy, delta] = self.findIntegerCoordinates(below(1), below(2), altitude, false, delta);
                below = [ix iy];
            else
                below = [];
            end
            if ~isempty(args.above)
                above = args.above;
                [ix, iy, delta] = self.findIntegerCoordinates(above(1), above(2), altitude, true, delta);
                above = [ix iy];
            else
                above = [];
            end
            [lip, delta] = self.lipFromBisection(altitude, below, above, delta);
        end

    end


    methods (Access = protected) % Utility methods

        function delta = powerOfTwoDelta(self, realDelta)
        % Returns the largest power-of-two integer delta smaller than the given real delta
        %
        % Raises:
        %   Raises an error if the corresponding delta is too small.
        %
        % Args:
        %   realDelta (double): Delta in real coordinates
        %
        % Returns:
        %   integer: Power-of-two delta in integer coordinates
            deltax = initialStepSize/abs(self.data.xRange(2) - self.data.xRange(1))*self.data.xDivisions;
            deltay = initialStepSize/abs(self.data.yRange(2) - self.data.yRange(1))*self.data.yDivisions;
            delta = 2^floor(log2(min(deltax, deltay)));
            assert(delta >= 1, 'Delta too small, augment xDivisions and yDivisions');
        end

    end

    methods (Access = protected) % Isoline methods

        function plotIsoline(self, j)
        % Updates the plot of the given isoline
            if ~isempty(self.figure) && self.plotIsolines
                [x, y] = self.data.isolinePath(self.evaluationFunction, j);
                figure(self.figure);
                h = self.figureIsolines{j};
                if isempty(h)
                    h = plot(x, y, 'b-');
                    self.figureIsolines{j} = h;
                else
                    set(h, 'XData', x, 'YData', y);
                end
                drawnow;
            end
        end

        function j = newIsoline(self, altitude, lip)
        % Creates a new isoline containing the given lip
        %
        % Args:
        %   f (function_handle): 2D function to plot
        %   altitude (double): Altitude of the new isoline
        %   lip (integer(1,4)): Integer coordinates describing a linearly interpolated point
        %
        % Returns:
        %   j: Index of the new isoline
            [isValid, reason] = self.data.isValidLip(f, altitude, lip);
            if ~isValid
                error(['Invalid lip: ' reason]);
            end
            j = length(self.data.altitudes + 1);
            self.data.altitudes(j) = altitude;
            self.data.paths{j} = lip(:);
            self.figureIsolines{j} = [];
            self.plotIsoline(j);
        end

    end

    methods (Access = protected) % Low-level methods

        function [ix, iy, delta] = findIntegerCoordinates(self, x, y, altitude, isAbove, delta0)
        % Finds a grid point below or above the given altitude near the given real coordinates
        %
        % Args:
        %   x (double): X real coordinate
        %   y (double): Y real coordinate
        %   altitude (double): Isoline altitude
        %   isAbove (logical): If true, finds a point where ``f(x,y) >= altitude``, else where ``f(x,y) < altitude``
        %   delta0 (integer): Starting power-of-two delta
        %
        % Returns
        % -------
        %   ix: integer
        %     Integer X coordinate
        %   iy: integer
        %     Integer Y coordinate
        %   delta: integer
        %     Updated delta
            f = self.evaluationFunction;
            delta = delta0;
            [ix0, iy0] = self.data.integerCoordinates(x, y);
            while delta >= 1
                for dx = 0:1
                    for dy = 0:1
                        ix = floor(ix0/delta+dx)*delta;
                        iy = floor(iy0/delta+dy)*delta;
                        c = self.data.eval(f, ix, iy) >= altitude;
                        if c == isAbove
                            return
                        end
                    end
                end
                delta = delta / 2;
            end
            error('Could not find starting point');
        end

    end

    methods (Access = protected) % Lip methods

        function [lip, delta] = lipFromBisection(self, altitude, below, above, delta0)
        % Finds a marching squares lip for the given altitude by grid bisection
        %
        % Args:
        %   altitude (double): Isoline altitude
        %   below (integer(1,2) or ``[]``): Known point below the altitude, optional
        %   above (integer(1,2) or ``[]``): Known point above the altitude, optional
        %   delta0 (integer): Starting power-of-two delta
        %
        % Returns
        % -------
        %   lip: integer(1,4)
        %     Linearly interpolated point
        %   delta: integer
        %     Updated delta
            f = self.evaluationFunction;
            squares = [0;0;self.data.xDivisions;self.data.yDivisions];
            i = 1;
            while i <= size(squares, 2)
                square = squares(:,i);
                x1 = square(1);
                y1 = square(2);
                x2 = square(3);
                y2 = square(4);
                pts = [x1 y1
                       x1 y2
                       x2 y1
                       x2 y2];
                for j = 1:4
                    if ~isempty(below) && ~isempty(above)
                        delta = min(delta0, gcd(gcd(below(1), below(2)), gcd(above(1), above(2))));
                        lip = self.lipFromArbitraryIntegerCoordinates(below, above, altitude, delta)
                        return
                    end
                    x = pts(j, 1);
                    y = pts(j, 2);
                    if self.data.eval(f, x, y) >= altitude
                        % above
                        if isempty(above)
                            above = [x y];
                        end
                    else
                        if isempty(below)
                            below = [x y];
                        end
                    end
                end
                if abs(x2-x1) > 1 && abs(y2-y1) > 1
                    x = (x1 + x2)/2;
                    y = (y1 + y2)/2;
                    assert(round(x) == x && round(y) == y);
                    newSquares = [x1 x  x1 x
                                  y1 y1 y  y
                                  x  x2 x  x2
                                  y  y  y2 y2];
                    squares = [squares newSquares];
                end
                i = i + 1;
            end
            error('Could not find a lip from grid bisection');
        end

        function lip = lipFromArbitraryIntegerCoordinates(self, below, above, altitude, delta)
        % Finds a marching squares lip for the given altitude, with integer points below and above as hints
        %
        % Args:
        %   below (integer(1,2)): Integer coordinates  of a point below the given altitude
        %   above (integer(1,2)): Integer coordinates of a point above the given altitude
        %   altitude (double): Isoline altitude
        %   delta (integer): Power-of-two delta
        %
        % Returns:
        %   lip: Linearly interpolated point
            f = self.evaluationFunction;
            bx = below(1);
            by = below(2);
            ax = above(1);
            ay = above(2);
            while (abs(ax - bx) > delta) || (abs(ay - by) > delta)
                x = floor((ax + bx)/2/delta)*delta;
                y = floor((ay + by)/2/delta)*delta;
                if self.data.eval(f, x, y) >= altitude
                    % the point is above
                    ax = x;
                    ay = y;
                else
                    bx = x;
                    by = y;
                end
            end
            assert((ax ~= bx) || (ay ~= by));
            if (ax == bx) || (ay == by)
                % already a proper lip
            else
                if self.data.eval(f, ax, by) >= altitude
                    ay = by;
                else
                    bx = ax;
                end
            end
            lip = [bx by ax ay];
        end

        function [lip, delta] = lipFromArbitraryRealCoordinates(self, below, above, altitude, delta0)
        % Finds a marching squares lip for the given altitude, with points below and above as hints
        %
        % Args:
        %   below (double(1,2)): Real coordinates of a point below the given altitude
        %   above (double(1,2)): Real coordinates of a point above the given altitude
        %   altitude (double): Isoline altitude
        %   delta0 (integer): Starting power-of-two delta
        %
        % Returns
        % -------
        %   lip: integer(1,4)
        %     Linearly interpolated point
        %   delta: integer
        %     Updated delta
            f = self.evaluationFunction;
            [bx, by, delta] = self.findIntegerCoordinates(below(1), below(2), altitude, false, delta0);
            [ax, ay, delta] = self.findIntegerCoordinates(above(1), above(2), altitude, true, delta);
            lip = self.lipFromArbitraryIntegerCoordinates([bx by], [ax ay], altitude);
        end

        function lip1 = splitLip(self, lip, altitude)
        % Splits a lip in half
        %
        % The lip must be parallel to either of the axes and have a delta >1 and a power-of-two.
        %
        % Args:
        %   lip (integer(1,4)): Integer coordinates describing a linearly interpolated point
        %   altitude (double): Altitude of the new isoline
        %
        % Returns:
        %   lip1: integer(1,4)
        %     Lip with half the delta
            f = self.evaluationFunction;
            ix1 = lip(1);
            iy1 = lip(2);
            ix2 = lip(3);
            iy2 = lip(4);
            delta = self.lipDelta(self, lip);
            assert(round(log2(delta)) == log2(delta), 'Delta must be a power-of-two');
            ix = (ix1 + ix2)/2;
            iy = (iy1 + iy2)/2;
            c1 = f(ix1, iy1) >= altitude;
            c2 = f(ix2, iy2) >= altitude;
            c = f(ix, iy) >= altitude;
            if c == c1
                ix1 = ix;
                iy1 = iy;
            else
                assert(c == c2);
                ix2 = ix;
                iy2 = iy;
            end
            lip1 = [ix1 iy1 ix2 iy2];
        end

    end

% $$$
% $$$         function splitLastLip(self, f)
% $$$         % Splits the last lip in half
% $$$         %
% $$$         % This method is used to resolve ambiguous situations such as when the four corners have the signs::
% $$$         %   +-  -+
% $$$         %   -+  +-
% $$$         %
% $$$         % This method modifies the last lip in place.
% $$$         %
% $$$         % Args:
% $$$         %   f (function_handle): Function ``f(x,y)`` to plot
% $$$             ox = self.path(1, end);
% $$$             oy = self.path(2, end);
% $$$             ix = self.path(3, end);
% $$$             iy = self.path(4, end);
% $$$             ic = self.eval(f, ix, iy) >= 0;
% $$$             oc = self.eval(f, ox, oy) >= 0;
% $$$             if ix == ox % if it is vertical lip
% $$$                 delta = oy - iy; % recover the step delta
% $$$                 assert(abs(delta) > 1, 'Error: grid too small, increase yDivisions and recompute.');
% $$$                 y = iy + delta/2;
% $$$                 c = self.eval(f, ix, y) >= 0;
% $$$                 if c == ic
% $$$                     % the new point is inside
% $$$                     iy = y;
% $$$                 else
% $$$                     % the new point is outside
% $$$                     oy = y;
% $$$                 end
% $$$             else
% $$$                 assert(iy == oy); % horizontal lip
% $$$                 delta = ox - ix; % recover the step delta
% $$$                 assert(abs(delta) > 1, 'Error: grid too small, increase yDivisions and recompute.');
% $$$                 x = ix + delta/2;
% $$$                 c = self.eval(f, x, iy) >= 0;
% $$$                 if c == ic
% $$$                     % the new point is inside
% $$$                     ix = x;
% $$$                 else
% $$$                     % the new point is outside
% $$$                     ox = x;
% $$$                 end
% $$$             end
% $$$             self.path(:, end) = [ox;oy;ix;iy];
% $$$         end
% $$$
% $$$         function step(self, f, dx, dy)
% $$$         % Performs a step of the boundary discovery
% $$$         %
% $$$         % When calling this function, a marching square is described implicitly by the last lip, which provides a
% $$$         % side of the square. The ``(dx, dy)`` argument is added to the two points of the lip, and describes the other
% $$$         % two corners of the square.
% $$$         %
% $$$         % Adds a lip to `.path`
% $$$         %
% $$$         % Args:
% $$$         %   f (function_handle): Function ``f(x,y)`` to plot
% $$$         %   dx (integer): X delta to add to the last lip (integer)
% $$$         %   dy (integer): Y delta to add to the last lip (integer)
% $$$             ox = self.path(1, end); % we extract the last lip
% $$$             oy = self.path(2, end);
% $$$             ix = self.path(3, end);
% $$$             iy = self.path(4, end);
% $$$             % we find the other corners of the marching square, the sides of the squares are the segments:
% $$$             % i-o (already in the path), o-a, i-b, a-b
% $$$             ax = ox + dx;
% $$$             ay = oy + dy;
% $$$             bx = ix + dx;
% $$$             by = iy + dy;
% $$$             % we compute/retrieve the signs of the four corners
% $$$             o = self.eval(f, ox, oy) >= 0;
% $$$             i = self.eval(f, ix, iy) >= 0;
% $$$             a = self.eval(f, ax, ay) >= 0;
% $$$             b = self.eval(f, bx, by) >= 0;
% $$$             % now we try all possible lips where the function has different signs
% $$$             n = 0;
% $$$             if o ~= a
% $$$                 n = n + 1;
% $$$                 ox1 = ox;
% $$$                 oy1 = oy;
% $$$                 ix1 = ax;
% $$$                 iy1 = ay;
% $$$             end
% $$$             if i ~= b
% $$$                 n = n + 1;
% $$$                 ox1 = bx;
% $$$                 oy1 = by;
% $$$                 ix1 = ix;
% $$$                 iy1 = iy;
% $$$             end
% $$$             if a ~= b
% $$$                 n = n + 1;
% $$$                 if a == o
% $$$                     ox1 = ax;
% $$$                     oy1 = ay;
% $$$                     ix1 = bx;
% $$$                     iy1 = by;
% $$$                 else
% $$$                     ox1 = bx;
% $$$                     oy1 = by;
% $$$                     ix1 = ax;
% $$$                     iy1 = ay;
% $$$                 end
% $$$             end
% $$$             if n > 1 % several solutions
% $$$                 warning('Ambiguous square, halving step size');
% $$$                 self.splitLastLip(f);
% $$$                 self.step(f, dx/2, dy/2);
% $$$             else % everything went well, add point
% $$$                 self.path(:,end+1) = [ox1;oy1;ix1;iy1];
% $$$             end
% $$$         end
% $$$
% $$$         function closePath(self, f)
% $$$         % Runs the marching squares algorithm
% $$$         %
% $$$         % This method assumes that `.initializePath` has already been called.
% $$$         %
% $$$         % It will find a path around the boundary to plot by applying the marching squares algorithm.
% $$$         %
% $$$         % Args:
% $$$         %   f (function_handle): Function ``f(x,y)`` to plot
% $$$             while any(self.path(:,1) ~= self.path(:,end))
% $$$                 % perform one step
% $$$                 ox = self.path(1, end);
% $$$                 oy = self.path(2, end);
% $$$                 ix = self.path(3, end);
% $$$                 iy = self.path(4, end);
% $$$                 if ix == ox % if it is vertical lip
% $$$                     delta = abs(oy - iy); % recover the step delta
% $$$                     prev1 = self.path(1, end-1);
% $$$                     prev2 = self.path(3, end-1);
% $$$                     % now we identify whether we came to the present lip from the left or the right
% $$$                     if prev1 < ix || prev2 < ix
% $$$                         % we came from the left, move to the right
% $$$                         self.step(f, delta, 0);
% $$$                     else
% $$$                         % we came from the right, move to the left
% $$$                         self.step(f, -delta, 0);
% $$$                     end
% $$$                 else
% $$$                     assert(iy == oy); % horizontal lip
% $$$                     delta = abs(ox - ix); % recover the step delta
% $$$                     prev1 = self.path(2, end-1);
% $$$                     prev2 = self.path(4, end-1);
% $$$                     % now we identify whether we came to the present lip from above or below
% $$$                     if prev1 < iy || prev2 < iy
% $$$                         % we came from above, move below
% $$$                         self.step(f, 0, delta);
% $$$                     else
% $$$                         % we came from below, move above
% $$$                         self.step(f, 0, -delta);
% $$$                     end
% $$$                 end
% $$$             end
% $$$         end
% $$$
% $$$     end

% $$$     methods % Plot retrieval
% $$$
% $$$
% $$$         function plotCachedPoints(self)
% $$$         % For debugging purposes, plots the points already computed and cached
% $$$         %
% $$$         % We plot the negative points in blue, and the nonnegative points in red
% $$$             [ix, iy, c] = find(self.data);
% $$$             [x, y] = self.realCoordinates(ix, iy);
% $$$             mask = c >= 0;
% $$$             plot(x(mask), y(mask), 'bx');
% $$$             plot(x(~mask), y(~mask), 'rx');
% $$$         end
% $$$
% $$$     end

% $$$         function delta = lipDelta(self, lip)
% $$$         % Returns the distance between the two points, given a lip parallel to either axis
% $$$         %
% $$$         % Args:
% $$$         %   lip (integer(1,4)): Integer coordinates describing a linearly interpolated point
% $$$         %
% $$$         % Returns:
% $$$         %   integer: Distance between the two points of the lip
% $$$             if lip(1) == lip(3) % vertical lip
% $$$                 delta = abs(lip(2) - lip(4));
% $$$             elseif lip(2) == lip(4) % horizontal lip
% $$$                 delta = abs(lip(1) - lip(3));
% $$$             else
% $$$                 error('Not a lip parallel to either the x or the y axis');
% $$$             end
% $$$         end
% $$$

    end

end
