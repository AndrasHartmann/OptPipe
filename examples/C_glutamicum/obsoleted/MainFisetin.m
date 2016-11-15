%Script for pelargonidin


%% Init
clc;
clear;
path=mfilename('fullpath');
[pathstr,name] = fileparts(path);
a=strfind(pathstr,'/dev');
%addpath(genpath(pathstr(1:a)));
%savepath;
a=strfind(pathstr,'/');
result_dir=strcat(pathstr(1:a(end)),'Results_glutamicum/fisetin');
mkdir(result_dir);


changeCobraSolver('gurobi5');

BM_rxn='biomass_a';
target_rxn='sink_fisetin';
max_KOs=03;


%% Read the model, insert heterologous reactions and constraints

model = readCbModel(strcat(pathstr,'/Glutamicum_models/ReliableModel_ElisabethZelle/iEZ475.xml'));
cobra_model=Get_heterologous(model,'fisetin');

%insert constraints from Zelle + coumaric acid (5 mM are added to the
%medium, assume import rate of 0.5)
cobra_model=Insert_constraints(cobra_model);
cobra_model=changeRxnBounds(cobra_model,'sink_coumarate',-0.5,'l');

cobra_model=changeObjective(cobra_model,target_rxn);
FBA.target=optimizeCbModel(cobra_model);

cobra_model=changeObjective(cobra_model,BM_rxn);
FBA.growth=optimizeCbModel(cobra_model);
maxBM=FBA.growth.f;


%% Get irrelevant reactions
irrelevantRxns(model)
load('irrelevantRxns.mat');
noDel=irrel.ind;
selectedRxnIDs=[1:1:length(model.rxns)];
selectedRxnIDs(noDel)=[];

clc;
%% OptKnock
disp('---------Optknock----------');

[optknock_results]=runOptknock(cobra_model,selectedRxnIDs,BM_rxn,target_rxn, max_KOs);

if sum(sum(~cellfun(@isempty,optknock_results)))==0
    disp('No KOs predicted by optknock')
end
fprintf('\n')
%% Robokod
disp('---------Robokod----------');

robokod_results=robokod(cobra_model, BM_rxn, target_rxn, max_KOs, selectedRxnIDs);
fprintf('\n')
%% OptGene
disp('---------OptGene----------');


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

disp('---------Enumeration with RobustKnock objective----------');
% Setting the objective for target
p = zeros(size(cobra_model.c));
targetRxnInd = find(ismember(cobra_model.rxns,  target_rxn));
p(targetRxnInd) = 1;
evalf = @(model) robustKnockSolution(model,p);

% Iterate through all combinations of deletions
[min_screening_results,fluxSolution_min] = testRxnDeletion(cobra_model, ... 
    evalf, max_KOs, selectedRxnIDs);


%Filter results
I=find(fluxSolution_min(:,1)>0.5*maxBM);
screening_i_min=min_screening_results(I,:);
screening_results_min=cobra_model.rxns(screening_i_min);


I=find(fluxSolution_max(:,1)>0.5*maxBM);
screening_i_max=max_screening_results(I,:);
screening_results_max=cobra_model.rxns(screening_i_max);

fprintf('\n')

%% Organize and analyze results

%optknock and optgene results may have only 2 columns
if size(optgene_results,2)<max_KOs
    optgene_results(:,end+1:max_KOs)={''};
end
if size(optknock_results,2)<max_KOs
    optknock_results(:,end+1:max_KOs)={''};
end

%assign method name to each result set
optknock_results(:,end+1)={'optknock'};
robokod_results(:,end+1)={'robokod'};
optgene_results(:,end+1)={'optgene'};
screening_results_max(:,end+1)={'screening max'};
screening_results_min(:,end+1)={'screening min'};

%join
%results is for output, all_comb is for use in the next step
results=vertcat(optknock_results,robokod_results,optgene_results,screening_results_max,screening_results_min);
all_comb=results(:,1:max_KOs);

%remove empty rows
a=sum(cellfun(@isempty,all_comb),2);
all_comb(a==3,:)=[];
results(a==3,:)=[];

%remove duplicate rows and register duplicates
[~,idx]=unique(strcat(all_comb(:,1),all_comb(:,2),all_comb(:,3)) , 'rows');
duplicate_i=setdiff(1:length(all_comb),idx);
duplicates=all_comb(duplicate_i,:);
all_comb=all_comb(idx,:);

%remove this part latter
all_comb(66,:)=[];
all_comb(201,:)=[];

%Get list of all reactions and analyze (subsystem enrichment)
[table_subsystems,all_rxns] = analyze_results(all_comb,model);

all_rxns(:,2)=printRxnFormula(model,all_rxns(:,1),false);



%% Rank combinations

disp('---------Ranking the KO combinations----------');

[FVA.min,FVA.max] = fluxVariability(cobra_model);

for i=1:length(results)
    comb=results(i,1:max_KOs);
    comb(~cellfun('isempty',comb));
    [biomass(i,1) target_min(i,1) target_max(i,1) distance(i,1)] = validateDeletion(cobra_model,...
        comb,target_rxn,FVA);
end


%filter out mutants that have max target=0 or biomass<0.5*maxBM
remove_ind=find(target_max==0);  
remove_ind=[remove_ind; find(biomass<0.5*maxBM)];

biomass(remove_ind)=[];
target_min(remove_ind)=[];
target_max(remove_ind)=[];
distance(remove_ind)=[];
results(remove_ind,:)=[]; 

%choose only mutants that have minimal target >0
choose_ind=find(target_min>0); 

best_results=results(choose_ind,:);
best_biomass=biomass(choose_ind);
best_target_min=target_min(choose_ind);
best_target_max=target_max(choose_ind);
best_distance=distance(choose_ind);

%organize
n=length(results);
table_results={'','','','method','biomass','minimal target','maximal target','distance to Wt'};
table_results(2:n+1,1:max_KOs+1)=results;
table_results(2:n+1,max_KOs+2)=num2cell(biomass);
table_results(2:n+1,max_KOs+3)=num2cell(target_min);
table_results(2:n+1,max_KOs+4)=num2cell(target_max);
table_results(2:n+1,max_KOs+5)=num2cell(distance);

n=length(best_results);
table_best_results={'','','','method','biomass','minimal target','maximal target','distance to Wt'};
table_best_results(2:n+1,1:max_KOs+1)=best_results;
table_best_results(2:n+1,max_KOs+2)=num2cell(best_biomass);
table_best_results(2:n+1,max_KOs+3)=num2cell(best_target_min);
table_best_results(2:n+1,max_KOs+4)=num2cell(best_target_max);
table_best_results(2:n+1,max_KOs+5)=num2cell(best_distance);

fprintf('\n')
%% Best candidates

disp('---------best candidates----------');
fprintf('\n')

[~,I] = sort(best_distance,'ascend');
a=best_results(I,:);
disp('Best candidates, taking into account adaptability:')
disp(a(1:5,:))

fprintf('\n')
[~,I] = sort(best_target_max,'descend');
a=best_results(I,:);
disp('Best candidates, taking into account maximal production:')
disp(a(1:5,:))

fprintf('\n')
[~,I] = sort(best_target_min,'descend');
a=best_results(I,:);
disp('Best candidates, taking into account minimal production:')
disp(a(1:5,:))


fprintf('\n')

%% Output results
disp('---------Exporting to excel----------');

cd(result_dir);


%For windows:
%xlswrite('fisetin_results',all_rxns,'All KOs') 
%xlswrite('fisetin_results.xls',table_subsystems,'Subsystems enrichment');
%xlswrite('fisetin_results.xls',table_all_comb,'All KO combinations');
%xlswrite('fisetin_results.xls',table_best_comb,'Best KO combinations');

%For mac:
%uses an external function from file exchange:
%http://www.mathworks.com/matlabcentral/fileexchange/37560-xlwrite---export-data-to-excel-from-matlab-on-mac-win

%This one works for linux as well
%http://www.mathworks.com/matlabcentral/fileexchange/37560-xlwrite---export-data-to-excel-from-matlab-on-mac-win



xlwrite('fisetin_results.xls',all_rxns,'All KO reactions'); %put formulas
xlwrite('fisetin_results.xls',table_subsystems,'Subsystems enrichment');
xlwrite('fisetin_results.xls',table_results,'All KO combinations');
xlwrite('fisetin_results.xls',table_best_results,'Best KO combinations');

disp('Document written');




