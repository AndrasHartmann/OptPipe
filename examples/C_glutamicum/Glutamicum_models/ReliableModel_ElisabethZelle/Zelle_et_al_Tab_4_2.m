%clear all;clc;
%changeCobraSolver('glpk','LP');changeCobraSolver('glpk','MILP');clear ans;
changeCobraSolver('gurobi5','LP');changeCobraSolver('gurobi5','MILP');clear ans;
 
M=readCbModel('iEZ475');
%clc

%% +Glucose:exponentiell phase
% % WT 
Model=M;

Model=changeRxnBounds(Model,'EX_glc(e)',-4.98,'u');
Model=changeRxnBounds(Model,'EX_glc(e)',-5.796,'l');
Model=changeRxnBounds(Model,'EX_o2(e)',-4.38,'l');
Model=changeRxnBounds(Model,'EX_o2(e)',-3.9,'u');
Model=changeRxnBounds(Model,'AC_d',-(0.79/15),'l');
Model=changeRxnBounds(Model,'AC_d',-(0.69/15),'u');
Model=changeRxnBounds(Model,'SUC_d',(0.47/15),'l');
Model=changeRxnBounds(Model,'SUC_d',(1.01/15),'u');
Model=changeRxnBounds(Model,'EX_lac_L(e)',(18/25),'u');%
Model=changeRxnBounds(Model,'EX_lac_L(e)',(8/25),'l');%
Model=changeRxnBounds(Model,'ACC_glycogen_c',-15,'l');%
lsgWT=optimizeCbModel(Model);

%%
% % DF1F0
Model=M;

Model=changeRxnBounds(Model,'ATPase',0,'b');
Model=changeRxnBounds(Model,'EX_glc(e)',-15.36,'u');
Model=changeRxnBounds(Model,'EX_glc(e)',-24,'l');
Model=changeRxnBounds(Model,'EX_o2(e)',-12.48,'l');
Model=changeRxnBounds(Model,'EX_o2(e)',-9.72,'u');
Model=changeRxnBounds(Model,'PYR_d',-(6.36/11),'l');
Model=changeRxnBounds(Model,'PYR_d',-(4.04/11),'u');
Model=changeRxnBounds(Model,'AC_d',-(2.83/23),'l');
Model=changeRxnBounds(Model,'AC_d',-(1.49/23),'u');
Model=changeRxnBounds(Model,'SUC_d',(1.31/21),'l');
Model=changeRxnBounds(Model,'SUC_d',(2.83/21),'u');
Model=changeRxnBounds(Model,'EX_lac_L(e)',(21.8/35),'u');
Model=changeRxnBounds(Model,'EX_lac_L(e)',(13.2/35),'l');
Model=changeRxnBounds(Model,'ACC_glycogen_c',-7,'l');
lsgDM=optimizeCbModel(Model);

lsgWT.f
lsgDM.f
