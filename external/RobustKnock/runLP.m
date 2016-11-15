function res = runLP(C, A, b_L, b_U, lb, ub)
%function res = runLP(C, A, b_L, b_U, lb, ub)
% Function to run a Linear Program (LP)

LPproblem.c = C;
LPproblem.osense = -1; % max
LPproblem.A = sparse(A);
LPproblem.lb = lb; LPproblem.ub = ub;

LPproblem.b = b_U;
%Equality constraints
LPproblem.csense(1:size(A,1),1) = 'E';
%TODO: should be upper / lower bound

res = solveCobraLP(LPproblem);
