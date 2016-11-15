function res = runMILP(C, A, B, lb, ub, intVars)
%function res = runMILP(C, A, B, lb, ub, intVars)
%Function to run Mixed Integer Linear Program (MILP)

MILPproblem.c = C;
MILPproblem.osense = -1; % max
MILPproblem.A = sparse(A);
MILPproblem.b_L = [];
MILPproblem.b_U = B;
MILPproblem.b = B;
MILPproblem.lb = lb; MILPproblem.ub = ub;
MILPproblem.x0 = [];
MILPproblem.vartype = char(ones(1,length(C)).*double('C'));
MILPproblem.vartype(intVars) = 'I';
MILPproblem.csense = char(ones(1,length(B)).*double('L'));

%%Some gurobi MILP parameters

%[MILPproblem, solverParams] = setParams(MILPproblem, true, solverParams); 
%disp('Run COBRA MILP')

% Using max number of processors
solverParams.Threads = java.lang.Runtime.getRuntime().availableProcessors;
solverParams.outputflag = 1;
%solverParams.Heuristics= 0.2;
solverParams.Presolve= 1;
%solverParams.PreDepRow= 1;
%solverParams.ConcurrentMIP=1;
solverParams.PreSparsify=1;
%For OptKnock
%solverParams.MIPFocus=3;
solverParams.MIPFocus=2;
%Robustknock
%solverParams.MIPFocus=1;
%For optKnock
solverParams.FeasibilityTol=1e-9;
%RobustKnock
%solverParams.FeasibilityTol=1e-7;
%solverParams.Cuts=3;
%solverParams.BranchDir=-1;
%solverParams.Cutoff=1;
solverParams.PrePasses=10000;
solverParams.ImproveStartTime=300;
solverParams.TimeLimit=500;
%optKnock
%solverParams.ImproveStartGap=10;
%Robustknock (maybe)
%solverParams.ImpliedCuts=2;

res = solveCobraMILP(MILPproblem,solverParams);
