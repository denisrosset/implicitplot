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
                axis([self.data.xRange self.data.yRange]);
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
                drawnow;
            end
        end

    end


    methods % Isoline creation

        function indices = isolinesUsingGrid(self, altitude, cellSize)
        % Finds isolines at a given altitude using a discretization scheme
        %
        % Args:
        %   altitude (double): Altitude at which to compute the isolines
        %   cellSize (double): Size of the grid cells; should be roughly half the size of the smallest component to find
        %
        % Returns:
        %   integer(1,\*): Indices of the created isolines
            f = self.evaluationFunction;
            delta = self.powerOfTwoDelta(cellSize);
            dX = self.data.xDivisions/delta;
            dY = self.data.yDivisions/delta;
            below = zeros(dX+1,dY+1);
            for i = 0:dX
                for j = 0:dY
                    below(i+1,j+1) = -(self.data.eval(f, i*delta, j*delta) < altitude);
                end
            end
            % code: -1 if below, 0 if above but unrecognized, {1,2,3} if recognized component
            comp = 1;
            indices = [];
            while 1
                [x, y] = find(~below, 1);
                if isempty(x)
                    return
                end
                x = x - 1;
                y = y - 1;
                plotDone = false;
                toCheck = [x; y];
                i = 1;
                while i >= 1
                    pt = toCheck(:, i);
                    i = i - 1;
                    x = pt(1);
                    y = pt(2);
                    below(x+1, y+1) = comp;
                    candidates = [x-1 x+1 x   x
                                  y   y   y-1 y+1];
                    for j = 1:4
                        x1 = candidates(1, j);
                        y1 = candidates(2, j);
                        if x1 >= 0 && x1 <= dX && y1 >= 0 && y1 <= dY
                            if below(x1+1, y1+1) == 0
                                i = i + 1;
                                toCheck(:, i) = [x1; y1];
                            elseif below(x1+1, y1+1) == -1 && ~plotDone
                                lip = [x y x1 y1]*delta;
                                indices(1,end+1) = self.isolineFromLip(altitude, lip);
                                plotDone = true;
                            end
                        end
                    end
                end
                comp = comp + 1;
            end
        end

        function j = isolineFromLip(self, altitude, lip)
            j = self.newIsoline(altitude, lip);
            [dx, dy] = self.isolineDirection([], lip);
            self.walkIsoline(j, dx, dy);
            if strcmp(self.data.isolineStatus(j), 'wip')
                self.data.paths{j} = fliplr(self.data.paths{j});
                if size(self.data.paths{j}, 2) == 1
                    self.walkIsoline(j, -dx, -dy);
                else
                    previousLip = self.data.paths{j}(:,end-1).';
                    currentLip = self.data.paths{j}(:,end).';
                    [dx, dy] = self.isolineDirection(previousLip, currentLip);
                    self.walkIsoline(j, dx, dy);
                end
            end
            assert(~strcmp(self.data.isolineStatus(j), 'wip'));
        end

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
                delta = max(min(self.data.xDivisions, self.data.yDivisions)/64, 1);
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
            lip = self.lipFromBisection(altitude, below, above, delta);
            j = self.isolineFromLip(altitude, lip);
        end

        function [dx, dy] = isolineDirection(self, previousLip, currentLip)
            x1 = currentLip(1);
            y1 = currentLip(2);
            x2 = currentLip(3);
            y2 = currentLip(4);
            if x1 == x2 % vertical lip
                delta = abs(y1 - y2);
                if isempty(previousLip)
                    if x1 == 0
                        dx = delta;
                    else
                        dx = -delta;
                    end
                else
                    if any(previousLip([1 3]) < x1) % we came from the left, move to the right
                        dx = delta;
                    else % we came from the right, move to the left
                        dx = -delta;
                    end
                end
                dy = 0;
            else
                delta = abs(x1 - x2);
                dx = 0;
                if isempty(previousLip)
                    if y1 == 0
                        dy = delta;
                    else
                        dy = -delta;
                    end
                else
                    if any(previousLip([2 4]) < y1) % we came from above, move below
                        dy = delta;
                    else % we came from below, move above
                        dy = -delta;
                    end
                end
            end
        end

        function reachedBorder = walkIsoline(self, j, dx, dy)
        % Adds lips to the end of an isoline until it reaches the border or loops back on itself
        %
        % Args:
        %   j (integer): Isoline index
        %   dx (integer): X delta to push the last lip on the isoline, must not push it on border
        %   dy (integer): Y delta to push the last lip on the isoline
        %
        % Returns:
        %   logical: Whether the isoline reached the figure border
            startLip = self.data.paths{j}(:,1).';
            previousLip = self.data.paths{j}(:,end).';
            currentLip = self.stepIsoline(j, dx, dy);
            if isempty(currentLip)
                reachedBorder = true;
                return
            end
            while any(startLip ~= currentLip) && any(startLip ~= currentLip([3 4 1 2]))
                [dx, dy] = self.isolineDirection(previousLip, currentLip);
                previousLip = currentLip;
                currentLip = self.stepIsoline(j, dx, dy);
                if isempty(currentLip)
                    reachedBorder = true;
                    return
                end
            end
            self.data.paths{j}(:,end+1) = startLip;
            reachedBorder = false;
        end

        function newLip = stepIsoline(self, j, dx, dy)
        % Performs a step to push the end of an isoline
        %
        % When calling this function, a marching square is described implicitly by the last lip, which provides a
        % side of the square. The ``(dx, dy)`` argument is added to the two points of the lip, and describes the other
        % two corners of the square.
        %
        % Args:
        %   j (integer): Isoline index
        %   dx (integer): X delta to add to the last lip (integer)
        %   dy (integer): Y delta to add to the last lip (integer)
        %
        % Returns:
        %   integer(1,4) or ``[]``: New lip or ``[]`` if we are at the boundary
            f = self.evaluationFunction;
            altitude = self.data.altitudes(j);
            lip = self.data.paths{j}(:,end).';
            x1 = lip(1); % we extract the last lip
            y1 = lip(2);
            x2 = lip(3);
            y2 = lip(4);
            % we find the other corners of the marching square, the sides of the squares are the segments:
            % 1-3 3-4 2-4 (1-2 is the current side)
            x3 = x1 + dx;
            y3 = y1 + dy;
            x4 = x2 + dx;
            y4 = y2 + dy;
            if ~all(self.data.validIntegerCoordinates([x3 x4], [y3 y4]))
                newLip = [];
                return
            end
            valid = 0; % number of valid lips
            lip = [];
            c1 = self.data.eval(f, x1, y1) >= altitude;
            c2 = self.data.eval(f, x2, y2) >= altitude;
            c3 = self.data.eval(f, x3, y3) >= altitude;
            c4 = self.data.eval(f, x4, y4) >= altitude;
            % side 1-3
            if c1 ~= c3
                valid = valid + 1;
                lip = [x1 y1 x3 y3];
            end
            if c3 ~= c4
                valid = valid + 1;
                lip = [x3 y3 x4 y4];
            end
            if c2 ~= c4
                valid = valid + 1;
                lip = [x2 y2 x4 y4];
            end
            if valid > 1 % several solutions
                warning('Ambiguous square, halving step size');
                lip1 = self.splitLip(lip, self.data.altitudes(j));
                self.data.paths{j}(:,end) = lip1;
                dx1 = dx/2;
                dy1 = dy/2;
                assert(max(dx1, dy1) >= 1, 'Pathological boundary');
                newLip = self.stepIsoline(j, dx1, dy1);
            else % everything went well, add point
                newLip = lip;
                self.data.paths{j}(:,end+1) = lip;
                self.plotIsoline(j);
                drawnow;
            end
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
            deltax = realDelta/abs(self.data.xRange(2) - self.data.xRange(1))*self.data.xDivisions;
            deltay = realDelta/abs(self.data.yRange(2) - self.data.yRange(1))*self.data.yDivisions;
            delta = 2^floor(log2(min(deltax, deltay)));
            assert(delta >= 1, 'Delta too small, augment xDivisions and yDivisions');
        end

    end

    methods (Access = protected) % Isoline methods

        function plotIsoline(self, j)
        % Updates the plot of the given isoline
            if ~isempty(self.figure) && self.plotIsolines
                [x, y] = self.data.isolinePath(self.evaluationFunction, j);
                h = self.figureIsolines{j};
                if isempty(h)
                    figure(self.figure);
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
        %   altitude (double): Altitude of the new isoline
        %   lip (integer(1,4)): Integer coordinates describing a linearly interpolated point
        %
        % Returns:
        %   j: Index of the new isoline
            f = self.evaluationFunction;
            [isValid, reason] = self.data.isValidLip(f, altitude, lip);
            if ~isValid
                error(['Invalid lip: ' reason]);
            end
            j = length(self.data.altitudes) + 1;
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
                        lip = self.lipFromArbitraryIntegerCoordinates(below, above, altitude, delta);
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
            delta = max(abs(ix1 - ix2), abs(iy1 - iy2));
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

end
