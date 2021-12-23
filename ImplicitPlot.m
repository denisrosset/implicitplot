classdef ImplicitPlot < handle
% Describes the 2D plot of an implicit equation f(x,y) == 0
%
% This class is optimized for implicit plots where the function evaluation is expensive. The computed function values
% are cached in a sparse matrix, and we maintain a linear correspondance between the plot range and the integer
% coordinates indexing the sparse matrix.
%
% The plot range is given by the properties `.xRange` and `.yRange`, which contain each two double values, the start and
% the end of the interval for the respective axes.
%
% For the x axis, the integer coordinates are between ``0`` and ``xDivisions`` included, where:
%
% - ``0`` maps to ``xRange(1)``
% - ``xDivisions`` maps to ``xRange(2)``.
%
% For the y axis, the integer coordinates are between ``0`` and ``yDivisions`` included, where:
%
% - ``0`` maps to ``yRange(1)``
% - ``yDivisions`` maps to ``yRange(2)``.
%
% In the sparse matrix, we have the following possible values. Here we write the integer coordinates ``(ix,iy)``
% corresponding to the point ``(x,y)``.
%
% - ``data(ix,iy) == 0`` means that ``f(x,y)`` hasn't been computed.
% - ``data(ix,iy) == realmin`` means that ``f(x,y) == 0``.
% - ``data(ix,iy) == c`` means that ``f(x,y) == c``.
%
% In MATLAB, the matrix indexing is one-based, so accesses to ``data(ix,iy)`` are shifted by one.
%
% Note that the value ``f(x,y) == realmin`` is encoded as if ``f(x,y) == 0``. ``realmin`` is the smallest positive
% normalized floating point number, ``2.2251e-308``.
%
% The path we plot is described by linearly interpolated points (lips).

% A lip is described by a pair of points where the function takes different signs.
% Say the two points have coordinates ``(x1, y1)`` and ``(x2, y2)``, and ``sign(f(x1,y1)) ~= sign(f(x2,y2))``.
% The lip is given by ``(x, y) == (1-t)*(x1, y1) + t*(x2, y2)`` where ``0 <= t <= 1`` solves the equation
% ``(1-t)*f(x1,y1) + t*f(x2,y2) == 0``.
%
% We ask that ``f(x1, y1)`` has the same sign as ``f(xRange(1), yRange(1))``, i.e. is "outside".
%
% Usually, a lip is situated on the side of a marching square. Of course, the coordinates ``(x1, y1)`` and ``(x2, y2)``
% are not given in floating-point format, rather as integer coordinates ``(ix1, iy1)`` and ``(ix2, iy2)``.
%
% The path we plot is thus described by the property `.plot`, where each column stores the ``(ix1, iy1, ix2, iy2)`` of
% a lip.

    properties (SetAccess = protected)
        xRange % (double(1,2)): Range of x values
        yRange % (double(1,2)): Range of y values
        data % (sparse double matrix): Evaluated function values
        xDivisions % (integer): Power of two describing the number of divisions across the x axis
        yDivisions % (integer): Power of two describing the number of divisions across the y axis
        path % (integer(4,N)): Path being plotted
    end

    methods (Static)

        function ip = empty(xRange, yRange, xDivisions, yDivisions)
        % Creates an empty implicit plot
        %
        % Args:
        %   xRange (double(1,2)): X axis range ``[xmin, xmax]``
        %   yRange (double(1,2)): Y axis range ``[ymin, ymax]``
        %   xDivisions (integer): Number of divisions of the x axis, default 2^16
        %   yDivisions (integer): Number of divisions of the y axis, default 2^16
            if nargin < 4
                yDivisions = 2^16;
            end
            if nargin < 3
                xDivisions = 2^16;
            end
            data = sparse(xDivisions+1, yDivisions+1);
            ip = ImplicitPlot(xRange, yRange, data, xDivisions, yDivisions, zeros(4, 0));
        end

        function ip = load(filename)
        % Reads an implicit plot from a MAT file
        %
        % Args:
        %   filename (charstring): Name of the MAT file to read
        %
        % Returns:
        %   `.ImplicitPlot`: An implicit plot
            s = load(filename);
            ip = ImplicitPlot.fromStruct(s);
        end

        function ip = fromStruct(s)
        % Creates an implicit plot from the data held in a struct
        %
        % This function is used for compatiblity with Octave, and to create MAT files that are interoperable with Python
        % and possibly other implementations of this plot.
            assert(s.version == 1);
            assert(strcmp(s.type, 'ImplicitPlot'));
            ip = ImplicitPlot(s.xRange, s.yRange, s.data, s.xDivisions, s.yDivisions, s.path);
        end

    end

    methods % Data output methods

        function save(self, filename)
        % Saves an implicit plot to a MAT file
        %
        % Args:
        %   filename (charstring): Name of the MAT file to write
            s = self.toStruct;
            save(filename, s);
        end

        function s = toStruct(self)
        % Extract the data from this implicit plot into a struct
        %
        % Returns:
        %   struct: Struct containing the data from which this implicit plot can be reconstructed
            s = struct('version', {1}, 'type', {'ImplicitPlot'}, 'xRange', {self.xRange}, 'yRange', {self.yRange}, ...
                       'data', {self.data}, 'xDivisions', {self.xDivisions}, 'yDivisions', {self.yDivisions});
        end

    end

    methods

        function self = ImplicitPlot(xRange, yRange, data, xDivisions, yDivisions, path)
        % Constructs an implicit plot object
        %
        % For the description of arguments, see the description of class properties
            self.xRange = xRange;
            self.yRange = yRange;
            self.data = data;
            self.xDivisions = xDivisions;
            self.yDivisions = yDivisions;
            self.path = path;
        end

        function c = eval(self, f, ix, iy)
        % Evaluates the given function at the given grid point, with caching
        %
        % Args:
        %   ix (integer): Grid x index
        %   iy (integer): Grid y index
        %   f (function_handle): Function to evaluate in case the point is not available in the cache
        %
        % Returns:
        %   double: Function value at the given grid point
            c = full(self.data(ix+1, iy+1));
            if c == 0
                [x, y] = self.realCoordinates(ix, iy);
                c = f(x, y);
                if c == 0
                    c = realmin;
                end
                self.data(ix+1, iy+1) = c;
            end
            if c == realmin
                c = 0;
            end
        end

        function [ix, iy] = integerCoordinates(self, x, y)
        % Converts real coordinates to integer coordinates
        %
        % Args:
        %   x (double(1,n)): X coordinate with ``xRange(1) <= x(i) <= xRange(2)`` for all ``i``
        %   y (double(1,n)): Y coordinate with ``yRange(1) <= y(i) <= yRange(2)`` for all ``i``
        %
        % Returns
        % -------
        %   ix(1,n): integer
        %     X integer coordinates
        %   iy(1,n): integer
        %     Y integer coordinates
            xRange = self.xRange;
            yRange = self.yRange;
            ix = (x - xRange(1))/(xRange(2) - xRange(1))*self.xDivisions;
            iy = (y - yRange(1))/(yRange(2) - yRange(1))*self.yDivisions;
            ix(x == xRange(1)) = 0;
            ix(x == xRange(2)) = self.xDivisions;
            iy(y == yRange(1)) = 0;
            iy(y == yRange(2)) = self.yDivisions;
        end

        function [x, y] = realCoordinates(self, ix, iy)
        % Converts integer coordinates to real coordinates
        %
        % Args:
        %   x (double(1,n)): X coordinate with ``0 <= x(i) <= self.xDivisions`` for all ``i``
        %   y (double(1,n)): Y coordinate with ``0 <= y(i) <= self.yDivisions`` for all ``i``
        %
        % Returns
        % -------
        %   x(1,n): double
        %     X real coordinate
        %   y(1,n): double
        %     Y real coordinate
            xRange = self.xRange;
            yRange = self.yRange;
            x = xRange(1) + (xRange(2) - xRange(1)) * ix/self.xDivisions;
            y = yRange(1) + (yRange(2) - yRange(1)) * iy/self.yDivisions;
            x(ix == self.xDivisions) = xRange(2); % to avoid numerical errors
            y(iy == self.yDivisions) = yRange(2); % to avoid numerical errors
        end

        function initializePath(self, f, insideX, insideY, initialStepSize)
        % Initializes the plot by finding a pair of lips
        %
        % One needs to provide a point ``(insideX, insideY)`` so that ``f(insideX, insideY)`` has a different sign than
        % the value of ``f`` at the boundary of the plot.
        %
        % This method will overwrite the `.path` property with the start of a path.
        %
        % Args:
        %   f (function_handle): 2D function to plot
        %   insideX (double): X coordinate of a point "inside"
        %   insideY (double): Y coordinate of a point "inside"
            deltax = initialStepSize/abs(self.xRange(2) - self.xRange(1))*self.xDivisions;
            deltay = initialStepSize/abs(self.yRange(2) - self.yRange(1))*self.yDivisions;
            delta = 2^floor(log2(min(deltax, deltay)));
            % find coordinates of inside point
            [ix, iy] = self.integerCoordinates(insideX, insideY);
            ix = floor((ix-1)/delta)*delta;
            iy = floor((iy-1)/delta)*delta;
            ic = self.eval(f, ix, iy);
            % outside point
            oy = 0;
            oc = self.eval(f, ix, oy);
            assert((oc >= 0) ~= (ic >= 0));
            while iy - oy > delta
                my = (oy + iy)/2;
                my = round(my/delta)*delta;
                mc = self.eval(f, ix, my);
                if (mc >= 0) == (oc >= 0)
                    % the middle point is outside
                    oy = my;
                    oc = mc;
                else
                    iy = my;
                    ic = mc;
                end
            end
            assert(iy-oy == delta);
            self.path = [ix;oy;ix;iy];
            self.step(f, delta, 0);
        end

        function step(self, f, dx, dy)
        % Performs a step of the boundary discovery
        %
        % When calling this function, a marching square is described implicitly by the last lip, which provides a
        % side of the square. The ``(dx, dy)`` argument is added to the two points of the lip, and describes the other
        % two corners of the square.
        %
        % Adds a lip to `.path`
        %
        % Args:
        %   f (function_handle): Function ``f(x,y)`` to plot
        %   dx (double): X delta to add to the last lip
        %   dy (double): Y delta to add to the last lip
            ox = self.path(1, end); % we extract the last lip
            oy = self.path(2, end);
            ix = self.path(3, end);
            iy = self.path(4, end);
            % we find the other corners of the marching square, the sides of the squares are the segments:
            % i-o (already in the path), o-a, i-b, a-b
            ax = ox + dx;
            ay = oy + dy;
            bx = ix + dx;
            by = iy + dy;
            % we compute/retrieve the signs of the four corners
            o = self.eval(f, ox, oy) >= 0;
            i = self.eval(f, ix, iy) >= 0;
            a = self.eval(f, ax, ay) >= 0;
            b = self.eval(f, bx, by) >= 0;
            % now we try all possible lips where the function has different signs
            done = false;
            if o ~= a
                ox1 = ox;
                oy1 = oy;
                ix1 = ax;
                iy1 = ay;
                done = true;
            end
            if i ~= b
                assert(~done);
                ox1 = bx;
                oy1 = by;
                ix1 = ix;
                iy1 = iy;
                done = true;
            end
            if a ~= b
                assert(~done);
                if a == o
                    ox1 = ax;
                    oy1 = ay;
                    ix1 = bx;
                    iy1 = by;
                else
                    ox1 = bx;
                    oy1 = by;
                    ix1 = ax;
                    iy1 = ay;
                end
            end
            self.path(:,end+1) = [ox1;oy1;ix1;iy1];
        end

        function closePath(self, f)
        % Runs the marching squares algorithm
        %
        % This method assumes that `.initializePath` has already been called.
        %
        % It will find a path around the boundary to plot by applying the marching squares algorithm.
            while any(self.path(:,1) ~= self.path(:,end))
                % perform one step
                ox = self.path(1, end);
                oy = self.path(2, end);
                ix = self.path(3, end);
                iy = self.path(4, end);
                if ix == ox % if it is vertical lip
                    delta = abs(oy - iy); % recover the step delta
                    prev1 = self.path(1, end-1);
                    prev2 = self.path(3, end-1);
                    % now we identify whether we came to the present lip from the left or the right
                    if prev1 < ix || prev2 < ix
                        % we came from the left, move to the right
                        self.step(f, delta, 0);
                    else
                        % we came from the right, move to the left
                        self.step(f, -delta, 0);
                    end
                else
                    assert(iy == oy); % horizontal lip
                    delta = abs(ox - ix); % recover the step delta
                    prev1 = self.path(2, end-1);
                    prev2 = self.path(4, end-1);
                    % now we identify whether we came to the present lip from above or below
                    if prev1 < iy || prev2 < iy
                        % we came from above, move below
                        self.step(f, 0, delta);
                    else
                        % we came from below, move above
                        self.step(f, 0, -delta);
                    end
                end
            end
        end

        function [x, y] = computePath(self, f)
        % Retrieves the path to plot in real coordinates by performing linear interpolation
        %
        % Args:
        %   f (function_handle): Function to plot (normally, when this function is called, all req. values are cached)
        %
        % Returns
        % -------
        %   x: double(1, N)
        %    X coordinates of the path to plot
        %   y: double(1, N)
        %    Y coordinates of the path to plot
            n = size(self.path, 2);
            x = zeros(1, n);
            y = zeros(1, n);
            for i = 1:n
                ix1 = self.path(1, i);
                iy1 = self.path(2, i);
                ix2 = self.path(3, i);
                iy2 = self.path(4, i);
                c1 = self.eval(f, ix1, iy1);
                c2 = self.eval(f, ix2, iy2);
                assert(sign(c1) ~= sign(c2));
                % (1-t)*c1 + t*c2 = 0
                % t = -c1/(c2-c1)
                t = -c1/(c2-c1);
                [x1, y1] = self.realCoordinates(ix1, iy1);
                [x2, y2] = self.realCoordinates(ix2, iy2);
                x(i) = (1-t)*x1 + t*x2;
                y(i) = (1-t)*y1 + t*y2;
            end
        end

        function plotCachedPoints(self)
        % For debugging purposes, plots the points already computed and cached
        %
        % We plot the negative points in blue, and the nonnegative points in red
            [ix, iy, c] = find(self.data);
            [x, y] = self.realCoordinates(ix, iy);
            mask = c >= 0;
            plot(x(mask), y(mask), 'bx');
            plot(x(~mask), y(~mask), 'rx');
        end

    end

end
