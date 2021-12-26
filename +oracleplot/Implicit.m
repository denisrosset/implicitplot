classdef Implicit < handle
% Describes the implicit plot of isolines of a function (``f(x,y) == cte``)
%
% This class caches computed function values, and maintains a list of computed isolines.
%
% **Plot range and discretization**
%
% The computed function values are cached in a sparse matrix, and we maintain a linear correspondance between the
% plot range and the integer coordinates indexing the sparse matrix.
%
% The plot range is given by the properties `.xRange` and `.yRange`, which are row vectors of two double values,
% the start and the end of the interval for each of the axes.
%
% Internally, we discretize these ``x`` and ``y`` ranges linearly.
%
% For the x axis, the integer coordinates are between ``0`` and `.xDivisions` included, where:
%
% - ``0`` maps to ``xRange(1)``
% - ``xDivisions`` maps to ``xRange(2)``.
%
% For the y axis, the integer coordinates are between ``0`` and `.yDivisions` included, where:
%
% - ``0`` maps to ``yRange(1)``
% - ``yDivisions`` maps to ``yRange(2)``.
%
% **Sparse matrix data storage**
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
% **Isolines**
%
% We maintain a list of computed isolines.
%
% Each isoline is described by an altitude ``v`` and a series of 2D coordinates ``(x,y)`` that approximately satisfy
% ``f(x,y) == v``.
%
% Each isoline is in one of the following states:
%
% * closed: if the first and last points of the isoline are identical
% * open: if the first and last points of the isoline are on the plot boundary but not identical
% * wip: if the first and last points are not identical and at least one is not on the plot boundary
%
% Each isoline corresponds to an integer matrix with ``N`` rows and four columns, where ``N`` is the number of points.
% Each row in this matrix stores four integers which together describe a linearly interpolated point.
%
% **Lips: linearly interpolated points**
%
% The points describing an isoline are linearly interpolated points (lips).
%
% Each lip is described by a pair of coordinates where the function ``f`` takes different signs. By "sign" we mean
% whether ``f(x,y) >= v`` or ``f(x,y) < v``.
%
% Say the two points have coordinates ``(x1, y1)`` and ``(x2, y2)``, and ``sign(f(x1,y1)) ~= sign(f(x2,y2))``.
%
% The lip is given by ``(x, y) == (1-t)*(x1, y1) + t*(x2, y2)`` where ``0 <= t <= 1`` solves the equation
% ``(1-t)*f(x1,y1) + t*f(x2,y2) == 0``.
%
% The coordinates ``(x1, y1)`` and ``(x2, y2)`` are not given in floating-point format, rather as integer coordinates
% ``(ix1, iy1)`` and ``(ix2, iy2)``.
%
% Each lip is thus described by four integers ``(ix1, iy1, ix2, iy2)``.

    properties (SetAccess = protected) % Plot and cached evaluations
        xRange % (double(1,2)): Range of x values
        yRange % (double(1,2)): Range of y values
        xDivisions % (integer): Power of two describing the number of divisions across the x axis
        yDivisions % (integer): Power of two describing the number of divisions across the y axis
        data % (sparse double matrix): Evaluated function values
    end

    properties (SetAccess = protected)
        altitudes % (double(1,I)): Altitudes of the isolines
        paths % (cell(1,I) integer(Ni,4)): Path of the isolines
    end

    methods (Static)

        function ip = empty(xRange, yRange, varargin)
        % Creates an empty implicit plot
        %
        % Args:
        %   xRange (double(1,2)): X axis range ``[xmin, xmax]``
        %   yRange (double(1,2)): Y axis range ``[ymin, ymax]``
        %
        % Keyword Args:
        %   xDivisions (integer): Number of divisions of the x axis, default 2^16
        %   yDivisions (integer): Number of divisions of the y axis, default 2^16
            args = struct('xDivisions', 2^16, 'yDivisions', 2^16);
            args = oracleplot.populateStruct(args, varargin);
            data = sparse(args.xDivisions+1, args.yDivisions+1);
            ip = oracleplot.Implicit(xRange, yRange, data, args.xDivisions, args.yDivisions, zeros(1, 0), cell(1, 0));
        end

        function ip = load(filename)
        % Reads an implicit plot from a MAT file
        %
        % Args:
        %   filename (charstring): Name of the MAT file to read
        %
        % Returns:
        %   `.Implicit`: An implicit plot
            s = load(filename);
            ip = oracleplot.Implicit.fromStruct(s);
        end

        function ip = fromStruct(s)
        % Creates an implicit plot from the data held in a struct
        %
        % This function is used for compatiblity with Octave, and to create MAT files that are interoperable with Python
        % and possibly other implementations of this plot.
        %
        % Args:
        %   s (struct): Structure containing the data
        %
        % Returns:
        %   `.Implicit`: An implicit plot
            assert(s.version == 1);
            assert(strcmp(s.type, 'oracleplot.Implicit'));
            ip = oracleplot.Implicit(s.xRange, s.yRange, s.data, s.xDivisions, s.yDivisions, s.altitudes, s.paths);
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
            s = struct('version', {1}, 'type', {'oracleplot.Implicit'}, ...
                       'xRange', {self.xRange}, 'yRange', {self.yRange}, ...
                       'data', {self.data}, 'xDivisions', {self.xDivisions}, 'yDivisions', {self.yDivisions}, ...
                       'altitudes', {self.altitudes}, 'paths', {self.paths});
        end

    end

    methods % Constructor

        function self = Implicit(xRange, yRange, data, xDivisions, yDivisions, altitudes, paths)
        % Constructs an implicit plot object
        %
        % For the description of arguments, see the description of class properties
            self.xRange = xRange;
            self.yRange = yRange;
            self.data = data;
            self.xDivisions = xDivisions;
            self.yDivisions = yDivisions;
            self.paths = paths;
            self.paths = paths;
        end

    end

    methods % Low-level functions

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

        function [ix, iy] = unroundedIntegerCoordinates(self, x, y)
        % Converts real coordinates to integer coordinates
        %
        % Args:
        %   x (double(1,n)): X coordinate with ``xRange(1) <= x(i) <= xRange(2)`` for all ``i``
        %   y (double(1,n)): Y coordinate with ``yRange(1) <= y(i) <= yRange(2)`` for all ``i``
        %
        % Returns
        % -------
        %   ix: double(1,n)
        %     X integer coordinates (not rounded!)
        %   iy: double(1,n)
        %     Y integer coordinates (not rounded!)
            xRange = self.xRange;
            yRange = self.yRange;
            ix = (x - xRange(1))/(xRange(2) - xRange(1))*self.xDivisions;
            iy = (y - yRange(1))/(yRange(2) - yRange(1))*self.yDivisions;
            ix(x == xRange(1)) = 0;
            ix(x == xRange(2)) = self.xDivisions;
            iy(y == yRange(1)) = 0;
            iy(y == yRange(2)) = self.yDivisions;
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
        %   ix: double(1,n)
        %     X integer coordinates (rounded!)
        %   iy: double(1,n)
        %     Y integer coordinates (rounded!)
            [ix, iy] = self.unroundedIntegerCoordinates(x, y);
            ix = round(ix);
            iy = round(iy);
        end

        function [x, y] = realCoordinates(self, ix, iy)
        % Converts integer coordinates to real coordinates
        %
        % Args:
        %   ix (integer(1,n)): X coordinate with ``0 <= x(i) <= self.xDivisions`` for all ``i``
        %   iy (integer(1,n)): Y coordinate with ``0 <= y(i) <= self.yDivisions`` for all ``i``
        %
        % Returns
        % -------
        %   x: double(1,n)
        %     X real coordinate
        %   y: double(1,n)
        %     Y real coordinate
            xRange = self.xRange;
            yRange = self.yRange;
            x = xRange(1) + (xRange(2) - xRange(1)) * ix/self.xDivisions;
            y = yRange(1) + (yRange(2) - yRange(1)) * iy/self.yDivisions;
            x(ix == self.xDivisions) = xRange(2); % to avoid numerical errors
            y(iy == self.yDivisions) = yRange(2); % to avoid numerical errors
        end

        function l = validIntegerCoordinates(self, ix, iy)
        % Returns whether the given integer coordinates are valid
        %
        % Args:
        %   ix (integer(1,n)): X integer coordinates
        %   iy (integer(1,n)): Y integer coordinates
        %
        % Returns:
        %   logical(1,n): Whether each of the given coordinates is valid
            vx = (ix == round(ix)) & (ix >= 0) & (ix <= self.xDivisions);
            vy = (iy == round(iy)) & (iy >= 0) & (iy <= self.yDivisions);
        end

        function l = integerCoordinatesOnBoundary(self, ix, iy)
        % Returns whether the given integer coordinates are on the boundary
        %
        % Args:
        %   ix (integer(1,n)): X coordinate with ``0 <= x(i) <= self.xDivisions`` for all ``i``
        %   iy (integer(1,n)): Y coordinate with ``0 <= y(i) <= self.yDivisions`` for all ``i``
        %
        % Returns:
        %   logical(1,n): Whether each of the given coordinates is on the boundary
            onx = (ix == 0) | (ix == self.xDivisions);
            ony = (iy == 0) | (iy == self.yDivisions);
            l = onx & ony;
        end

    end

    methods % Linearly interpolated points manipulation

        function [x, y] = computeLipCoordinates(self, f, altitude, lip)
        % Computes the real coordinates corresponding to a linearly interpolated point
        %
        % Args:
        %   f (function_handle): 2D function to plot
        %   altitude (double): Altitude of the new isoline
        %   lip (integer(1,4)): Integer coordinates describing a linearly interpolated point
        %
        % Returns
        % -------
        %   x: double
        %     Interpolated X coordinate
        %   y: double
        %     Interpolated Y coordinate
            ix1 = lip(1);
            iy1 = lip(2);
            ix2 = lip(3);
            iy2 = lip(4);
            v1 = self.eval(f, ix1, iy1);
            v2 = self.eval(f, ix2, iy2);
            assert((v1 >= altitude) ~= (v2 >= altitude), 'Points of a lip must differ in position w.r.t. altitude');
            % We solve the equation
            % (1-t)*v1 + t*v2 = a
            % t*(v2-v1) + v1 = a
            % t*(v2-v1) = a-v1
            % t=(a-v1)/(v2-v1)
            t = (altitude-v1)/(v2-v1);
            [x1, y1] = self.realCoordinates(ix1, iy1);
            [x2, y2] = self.realCoordinates(ix2, iy2);
            x = (1-t)*x1 + t*x2;
            y = (1-t)*y1 + t*y2;
        end

        function [isValid, reason] = isValidLip(self, f, altitude, lip)
        % Returns whether the given lip is valid
        %
        % Args:
        %   f (function_handle): 2D function to plot
        %   altitude (double): Altitude of the new isoline
        %   lip (integer(1,4)): Integer coordinates describing a linearly interpolated point
        %
        % Returns
        % -------
        %   isValid: logical
        %     True if the lip is valid
        %   reason: charstring
        %     If the lip is valid, empty. If invalid, provides the reason.
            ix1 = lip(1);
            iy1 = lip(2);
            ix2 = lip(3);
            iy2 = lip(4);
            if any(~self.validIntegerCoordinates([ix1 ix2], [iy1 iy2]))
                isValid = false;
                reason = 'Point outside the plot range';
                return
            end
            c1 = self.eval(f, ix1, iy1) >= altitude;
            c2 = self.eval(f, ix2, iy2) >= altitude;
            if c1 == c2
                isValid = false;
                if c1
                    reason = 'Both points are above the isoline altitude';
                else
                    reason = 'Both points are below the isoline altitude';
                end
                return
            end
            isValid = true;
            reason = '';
        end

    end

    methods % Isoline manipulation

        function status = isolineStatus(self, f, j)
        % Returns the status of the given isoline
        %
        % Args:
        %   f (function_handle): 2D function to plot
        %   j (integer): Index of the isoline
        %
        % Returns:
        %   {'wip', 'open', 'closed'}: Isoline status
            if size(self.paths{j}, 1) == 1
                status = 'wip';
                return
            end
            first = self.paths{j}(1,:);
            last = self.paths{j}(end,:);
            if all(first == last)
                status = 'closed';
            else
                if all(self.integerCoordinatesOnBoundary([first last]))
                    status = 'open';
                else
                    status = 'wip';
                end
            end
        end

        function [x, y] = isolinePath(self, f, j)
        % Returns the real coordinates corresponding to the path of an isoline
        %
        % The returned data can be passed directly to MATLAB's ``plot`` function.
        %
        % Args:
        %   f (function_handle): 2D function to plot
        %   j (integer): Index of the isoline
        %
        % Returns
        % -------
        %   x: double(1, N)
        %    X coordinates of the path to plot
        %   y: double(1, N)
        %    Y coordinates of the path to plot
            path = self.paths{j};
            altitude = self.altitudes(j);
            n = size(path, 1);
            x = zeros(1, n);
            y = zeros(1, n);
            for i = 1:n
                lip = path(i, :);
                [xi, yi] = self.computeLipCoordinates(self, f, altitude, lip);
                x(i) = xi;
                y(i) = yi;
            end
        end

    end

end
