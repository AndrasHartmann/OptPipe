function [results idx]=runOptknock(model,selectedRxnIDs,BM_rxn,target_rxn,max_KO,minbiomass)

    %% Init
    optKnockFlag=1;
    robustKnockFlag=0;
    maxW=1000;

    objectiveInd = find(ismember(model.rxns, BM_rxn));
    chemicalInd = find(ismember(model.rxns, target_rxn));

    %%We only support gurobi
    changeCobraSolver('gurobi5','MILP');

    %% Exending model with some additional params
    [metNum, rxnNum] = size(model.S);
    model.row_lb=zeros(metNum,1);
    model.row_ub=zeros(metNum,1);

    model.C_chemical=zeros(rxnNum,1);
    model.C_chemical(chemicalInd)=1;

    model.organismObjective=zeros(rxnNum,1);
    model.organismObjective(objectiveInd)=1;

    model.organismObjectiveInd=objectiveInd;
    model.chemicalInd=chemicalInd;

    model.biomassInd=objectiveInd;

    %TODO: set biomass
    model.lb(objectiveInd)=minbiomass;
                                            
    res =  robustKnock(model, chemicalInd, objectiveInd, selectedRxnIDs, max_KO, maxW, optKnockFlag, robustKnockFlag);
    idx = res.Result_optKnock.knockedV';
    results = model.rxns(idx)';

end




