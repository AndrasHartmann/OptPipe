function success = Main()
%function Main()
% Main function for C. glutamicum optimizing naringenin production
% OUTPUT:
%        success of the identification using the pipeline

success = false;

branch = 'naringenin'; 

%% Init
clc;

%Disable warnings
warning('off', 'all');
%TODO: warnings should be handled

% Read the model, insert heterologous reactions and constraints
model = readCbModel([cd filesep 'Glutamicum_models' filesep 'ReliableModel_ElisabethZelle' filesep 'iEZ475.xml']);
[cobra_model tmp]=Get_heterologous(model, 'pelargonidin');
sink = 'sink_naringenin';

% parameter settings
BM_rxn='biomass_a';
target_rxn= sink;
max_KOs=2;

% currently only support gurobi solver
changeCobraSolver('gurobi5');

%insert constraints from Zelle + coumaric acid (5 mM are added to the
%medium, assume import rate of 0.5)
cobra_model=Insert_constraints(cobra_model);
cobra_model=changeRxnBounds(cobra_model,'sink_coumarate',-0.5,'l');

[optknock_results,robokod_results,optgene_results,screening_results_max,screening_results_min] = ...
optPipe(cobra_model,BM_rxn,target_rxn,max_KOs,branch)

success = true;

%Re-enable warnings
warning('on', 'all');
