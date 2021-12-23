function res = minCoeff(P)
% Returns the minimum coefficient of the given distribution
%
% Args:
%   P (double(\*,\*,...)): Probability distribution
%
% Returns:
%   double: Minimal coefficient in ``P``
    res = min(P(:));
end
