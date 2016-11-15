%adding path for common functions
addpath('../../common_functions');
%% Import model
out.model = readCbModel('../../Start_glutamicum/Glutamicum_models/ReliableModel_ElisabethZelle/iEZ475.xml');  
% add the sink
out.model = AddSinkToModel(out.model);
% add constraints
out.model = Insert_constraints(out.model)

%% Initialization
%Find biomass and McoA rxns
biomassrxn = find(ismember(out.model.rxns,'biomass_a')); %301
targetRxn = find(ismember(out.model.rxns,  'sink_malonylCoA')); % MCoA

% lower flux for reaction to be lethal 
fluxTol = 1e-9; 

% set cobra solver
% is this still needed?
%changeCobraSolver('gurobi6');

%% pre-processing: irrelevant reactions
[irrel relevantRxns] = irrelevantRxns(out.model);
irrelevant = irrel.blockedRxns|irrel.tr_synthRxns|irrel.lethalRxns;

%% Objective for target
p = zeros(size(out.model.c));
p(targetRxn) = -1;
%p(targetRxn) = 1;

evalf = @(model) robustKnockSolution(model,p);

%% Iterate through all combinations of deletions
[delRxn,fluxSolution] = testRxnDeletion(out.model, evalf, 3, ... 
    find(~irrelevant));

%% Display and plotting
for i = 1:size(delRxn, 1)
    disp('deletion');
    out.model.rxns(delRxn(i,:))'
    disp('prediction (biomass, target)');
    fluxSolution(i,:)
end

scatter(fluxSolution(:,1), fluxSolution(:,2))
xlabel('Biomass')
ylabel('Target')
