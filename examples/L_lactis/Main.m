function success = Main(branch)
%function Main(branch)
% Main function for L. lactis
% INPUT:
%       branch : string for identifying the branch, e.g. one of the four branches 
%                'resveratrol', 'pelargonidin', 'fisetin', 'quercetin'
% OUTPUT:
%        success of the identification using the pipeline


success = false;
if ~ismember(branch, { 'resveratrol', 'pelargonidin', 'fisetin', 'quercetin'})
    error('The parameter branch should be one of the four branches: ''resveratrol'', ''pelargonidin'', ''fisetin'', ''quercetin''');
    return;
end


%% Init
clc;

%Disable warnings
warning('off', 'all');
%TODO: warnings should be handled


% Read the model, including heterologous reactions and constraints
cobra_model = readCbModel([cd filesep 'lactis_models' filesep 'constrained_models' filesep 'lactis_' branch '.xml']);

BM_rxn='biomass_LLA';
%target_rxn='sink_dihydrofisetin';
%target_rxn='ACCOAC'
max_KOs=03;

% currently only support gurobi solver
changeCobraSolver('gurobi5');

switch branch
    case 'resveratrol' 
        target_rxn = 'sink_transresveratrol';
    case 'pelargonidin' 
        target_rxn = 'sink_pelargonidin';
    case 'fisetin' 
        target_rxn = 'sink_fisetin';
    case 'quercetin'
        target_rxn = 'sink_quercetin';
end;

%Set the objective (just in case)
cobra_model = changeObjective(cobra_model, BM_rxn)


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

[optknock_results,robokod_results,optgene_results,screening_results_max,screening_results_min] = ...
optPipe(cobra_model,BM_rxn,target_rxn,max_KOs,branch)

success = true;

%Re-enable warnings
warning('on', 'all');
