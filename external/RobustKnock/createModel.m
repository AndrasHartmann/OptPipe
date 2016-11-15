function consModel=createModel(model, chemicalInd, organismObjectiveInd, fixedGlu)

%lets start!
%load('Ecoli_iJR904.mat');
%consModel=Ecoli_iJR904;

%load('C_glutamicum');
%consModel = model_sink;
consModel = model;

[metNum, rxnNum] = size(consModel.S);
consModel.row_lb=zeros(metNum,1);
consModel.row_ub=zeros(metNum,1);
consModel.C_chemical=zeros(rxnNum,1);
consModel.C_chemical(chemicalInd)=1;
consModel.organismObjective=zeros(rxnNum,1);
consModel.organismObjective(organismObjectiveInd)=1;

consModel.organismObjectiveInd=organismObjectiveInd;
consModel.chemicalInd=chemicalInd;
medium = strmatch('EX_glc(e)',consModel.rxns);

%set medium constraints
consModel.lb( medium ) = -fixedGlu;
consModel.ub( medium ) = -fixedGlu;

biomassInd = [strmatch('Biomass', consModel.rxns) , strmatch('Ec_biomass', consModel.rxns), strmatch('biomass_a', consModel.rxns) ];
    if(length(biomassInd) < 1),
        disp('No Biomass reaction found!');
        return;
    end
    if(length(biomassInd) > 1),
        disp('Too many Biomass reactions found!');
        return;
    end
% set minimum growth constraint
consModel.biomassInd=biomassInd;
consModel.lb(biomassInd)=0.1;

%{
%remove small  reactions
ind=find(consModel.S>-(10^-3) & consModel.S<0);
ind2=find(consModel.S<(10^-3) & consModel.S>0);
consModel.S(ind)=0;
consModel.S(ind2)=0;
%}


