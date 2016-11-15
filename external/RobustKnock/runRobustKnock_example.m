% Script to run the original example from Tepper and Shlomi (2009)

%parameters
optKnockFlag=1;
robustKnockFlag=0;
knockoutNum=2;
maxWrobust=1e7;
%E coli
objectiveInd = 150; %biomass

%chemicalInd = 329; % ethanol
%chemicalInd = 409; % succinate
%chemicalInd = 368; % D-Lactate
%chemicalInd = 359; % H+
chemicalInd = 294; % Acetate
%chemicalInd = 331; % Formate
%chemicalInd = 335; % Fumerate

%C. glutamicum
%objectiveInd = 301; %biomass
%chemicalInd = 547; % MCoA

fixedGlu = 10;

%constants
maxWopt=1000;
if (robustKnockFlag)
    maxW=maxWrobust;
else
    maxW=maxWopt;
end

%maxW=1e7;
load('Ecoli_iJR904.mat');
model=Ecoli_iJR904;



consModel=createModel(model, chemicalInd, objectiveInd, fixedGlu);
reaction_indices = 1:length(consModel.rxns);
excluded_reactions = findExcRxns(consModel,1);
reaction_indices = setdiff(reaction_indices, [find(excluded_reactions)' chemicalInd objectiveInd]);

results =  robustKnock(consModel, chemicalInd, objectiveInd, reaction_indices, knockoutNum, maxW, optKnockFlag, robustKnockFlag);

