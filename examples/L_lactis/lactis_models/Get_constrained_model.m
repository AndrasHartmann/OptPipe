%% Init
clc;
clear;
changeCobraSolver('gurobi5');
model=readCbModel('253_2013_5140_MOESM6_ESM.xml');
new_model=model;

%% Disable model constraints
load('disable_constraints.mat');
for i=1:length(dis_constraints.val)
   name=dis_constraints.label{i,1};
   index=find(ismember(model.rxnNames,name));
   rxn=model.rxns{index,1};
   new_model=changeRxnBounds(new_model,{rxn,rxn},[dis_constraints.val(i,1),dis_constraints.val(i,2)],{'l','u'}); 
end


%% Insert constraints from supp file 5 (for growth rate=0.5)
load('constraints_05.mat'); 
for i=1:length(constraints.val)
   name=constraints.label{i,1};
   index=find(ismember(model.rxnNames,name));
   rxn=model.rxns{index,1};
   new_model=changeRxnBounds(new_model,{rxn,rxn},[constraints.val(i,1),constraints.val(i,2)],{'l','u'});   
end

%% Constrain ATPm (constant value!)
new_model=changeRxnBounds(new_model,'ATPM',0.92,'b');
new_model=changeObjective(new_model,'biomass_LLA');

%% Remove constraints for lactate and succinate
new_model=changeRxnBounds(new_model,{'EX_lac_L(e)','EX_lac_L(e)'},[0,1000],{'l','u'});
new_model=changeRxnBounds(new_model,'EX_succ(e)',1000,'u');

%% LDH mutant
new_model=changeRxnBounds(new_model,'LDH_L',0,'b');

%% Write
writeCbModel(new_model,'sbml','lactis_constrained.xml')