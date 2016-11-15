%THIS SCRIPT is obsolated, do not use

%Script for resveratrol in lactis

%% Init
clc;
clear;
changeCobraSolver('gurobi5');
[pathstr,name] = fileparts(which('MainLactis'));
a=strfind(pathstr,'/');
result_dir=strcat(pathstr(1:a(end)),'results/resveratrol_Lactis');
mkdir(result_dir);

BM_rxn='biomass_LLA';
target_rxn='sink_transresveratrol';
max_KOs=03;


%% Read model, insert heterolgous pathways and constraints
cd(pathstr);

%This model has already the constraints from Flahault,2013 for gr=0.5
%the constraints for lactate and succinate were removed
%already ldh mutant
%The script that was used to obtain this model is 'Get_constrained_model'
model=readCbModel('lactis_models/lactis_constrained.xml');


%Allow oxygen
% new_model=changeRxnBounds(new_model,'EX_o2(e)',-100,'l');
% new_model=changeRxnBounds(new_model,'O2t',100,'u');
%this part is left out: resveratrol can be produced in anaerobiosis


cobra_model=Get_heterologousLactis(model,'resveratrol');
cobra_model=changeObjective(cobra_model,'sink_transresveratrol');
fba=optimizeCbModel(cobra_model)

%% Get irrelevant reactions
irrelevantRxns_Lactis(model);
load('irrelevantRxns_Lactis.mat');
noDel=irrel.ind;
selectedRxnIDs=[1:1:length(cobra_model.rxns)];
selectedRxnIDs(noDel)=[];

clc;

%% OptKnock
disp('---------Optknock----------');

[optknock_results]=runOptknock(cobra_model,selectedRxnIDs,BM_rxn,target_rxn,max_KOs);

if sum(sum(~cellfun(@isempty,optknock_results)))==0
    disp('No KOs predicted by optknock')
end
fprintf('\n')

%% Robokod
disp('---------Robokod----------');

try
robokod_results=robokod(cobra_model, BM_rxn, target_rxn, max_KOs, selectedRxnIDs);
fprintf('\n')
catch me
    warning(['no Robokod results found: ' me.message]);
    robokod_results = [];
end

%ERROR WITH THE CONSTRAINTS

%% OptGene
disp('---------OptGene----------');


try
%Change model so that the "irrelevant reactions" don't have a gene
%association
out_model=cobra_model;
for i=1:length(selectedRxnIDs)
    in=selectedRxnIDs(i);
    if in <= length(model.grRules)
        out_model.grRules{in}='';
        out_model.rxnGeneMat(in,:)=0;
    end
end


%Export model
writeCbModel(cobra_model,'sbml',strcat(result_dir,'/','fisetin_model'));
disp('model created in results folder');


%Run optgene externally
disp('Run OptGene in OptFlux for gene deletion:');
disp('MOMA, 5000 eval, max 3 KOs, SPEA2, Gene Del, maximize BPCY with substrate as coumarate');


%Read optgene results
disp('Insert xlsx file with optgene results in the results directory') %'Example: 'rxn1|rxn2|rxn3');
%Ask for user input( if the results for optgene are already in the folder)
%to carry on
uiwait(msgbox('Is the optgene_<target>.xlsx file already in the results folder?'));
[num,optgene_results,raw] = xlsread(strcat(result_dir,'/optgene_fisetin.xlsx'));



%Transform genes in reactions
[l,w]=size(optgene_results);
[row,col]=find(model.rxnGeneMat);
for i=1:l
    for j=1:w
        gene=optgene_results{i,j};
        if ~isempty(gene)
            in=find(ismember(model.genes,gene));
            rxn=model.rxns{row(find(col==in)),1};
            optgene_results{i,j}=rxn;
        end
    end
end
%[results a]= findRxnsFromGenes(model, optgene_results(1,:),'listresults',true)

catch me
    warning(['no OptGene results found: ' me.message]);
    optgene_results = [];
end

%% Enumeration methods

%OptKnock target: What is the maximum Target when Maximizing Biomass?

uiwait(msgbox('Enumeration can take long time, make sure that you have activated parallel toolbox'));
disp('---------Enumeration with OptKnock objective----------');
%Ask for user input( if the results for optgene are already in the folder)
%to carry on

addpath('../methods/enumerate');

% Setting the objective for target
p = zeros(size(cobra_model.c));
targetRxnInd = find(ismember(cobra_model.rxns,  target_rxn));
p(targetRxnInd) = -1;
evalf = @(model) robustKnockSolution(model,p);

% Iterate through all combinations of deletions
[max_screening_results,fluxSolution_max] = testRxnDeletion(cobra_model, ... 
    evalf, max_KOs, selectedRxnIDs);


%RobustKnock target: What is the minimum Target when Maximizing Biomass?

