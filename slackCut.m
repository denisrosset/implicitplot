function slack = slackCut(P)
% Computes the slack for the cut inflation in the triangle
%
% The cut is done between A and C, so the inflation is A1-B1-C1
%
% Requires YALMIP.
%
% Args:
%   P (double(2,2,2)): Probability distribition
%
% Returns:
%   double: Negative values represent infeasibility
    P_ABC = reshape(P, [2 2 2]);
    P_AC = sum(P_ABC, 2);
    P_AB = sum(P_ABC, 3);
    P_BC = sum(P_ABC, 1);
    P_A = sum(P_ABC, [2 3]);
    P_C = sum(P_ABC, [1 2]);
    P_A_C = kron(P_C(:), P_A(:));
    P_A1B1C1 = sdpvar(8,1);
    P_A1B1C1 = reshape(P_A1B1C1, [2 2 2]);
    P_A1C1 = sum(P_A1B1C1, 2);
    P_A1B1 = sum(P_A1B1C1, 3);
    P_B1C1 = sum(P_A1B1C1, 1);
    slack = sdpvar;
    CONS = [P_A1C1(:) == P_A_C(:)
            P_A1B1(:) == P_AB(:)
            P_B1C1(:) == P_BC(:)
            P_A1B1C1(:) >= slack
            P_ABC(:) >= slack];
    optimize(CONS, -slack);
    slack = double(slack);
end
