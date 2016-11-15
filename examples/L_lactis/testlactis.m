%init
clear;

branch = 'resveratrol'

success = false;
if ~ismember(branch, { 'resveratrol', 'pelargonidin', 'fisetin', 'quercetin'})
    error('The parameter branch should be one of the four branches: ''resveratrol'', ''pelargonidin'', ''fisetin'', ''quercetin''');
    return;
end



%Disable warnings
%warning('off', 'all');
%TODO: warnings should be handled


% Read the model, including heterologous reactions and constraints
model = readCbModel([cd filesep 'lactis_models' filesep 'constrained_models' filesep 'lactis_' branch '.xml']);

biomassRxn='biomass_LLA';
%target_rxn='sink_dihydrofisetin';
%target_rxn='ACCOAC'
max_KOs=03;

model.c(1) = 1;

% currently only support gurobi solver
changeCobraSolver('gurobi5');

switch branch
    case 'resveratrol' 
        targetRxn = 'sink_transresveratrol';
    case 'pelargonidin' 
        targetRxn = 'sink_pelargonidin';
    case 'fisetin' 
        targetRxn = 'sink_fisetin';
    case 'quercetin'
        targetRxn = 'sink_quercetin';
end;

selectedRxnList = getSelectedRxns(model, biomassRxn)
selectedRxnIDs = find(ismember(model.rxns,selectedRxnList));

% Setting the objective for target
p = zeros(size(model.c));
targetRxnInd = find(ismember(model.rxns,  targetRxn));
p(targetRxnInd) = 1;
evalf = @(model) robustKnockSolution(model,p);

%Iterate through all combinations of deletions
[max_screening_results,fluxSolution_max] = testRxnDeletion(model, evalf, max_KOs, selectedRxnIDs);


%% Read model, insert heterolgous pathways and constraints
%cd(pathstr);

%This model has already the constraints from Flahault,2013 for gr=0.5
%the constraints for lactate and succinate were removed
%already ldh mutant
%The script that was used to obtain this model is 'Get_constrained_model'
%cobra_model=readCbModel('lactis_models/lactis_Fisetin.xml');

%cobra_model=Get_constrained_model(cobra_model)
%cobra_model=changeObjective(cobra_model,BM_rxn);


%Allow oxygen
% new_model=changeRxnBounds(new_model,'EX_o2(e)',-100,'l');
% new_model=changeRxnBounds(new_model,'O2t',100,'u');
%this part is left out: resveratrol can be produced in anaerobiosis


%cobra_model=Get_heterologousLactis(model,'resveratrol');
%cobra_model=changeObjective(cobra_model,'sink_transresveratrol');

%productionEnvelope(cobra_model,{}, 'b', target_rxn,BM_rxn,0,20)
