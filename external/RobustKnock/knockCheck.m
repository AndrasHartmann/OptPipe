function res=knockCheck(model, knockedInds, coupledFlag, coupled)
%if we use the seperated model, we need to knock the coupled reactions too
res=[];
if(coupledFlag)
    for k=1:length(coupled)
        if (find(knockedInds==coupled(k,1)))
            knockedInds(end+1)=coupled(k,2);
        end
    end
end

%find fluxes
ubNew=model.ub;
ubNew(knockedInds)=0;
lbNew=model.lb;
lbNew(knockedInds)=0;

%find optimal objective and then force it
%ProbMaxBiomass=lpAssign(-model.organismObjective, model.S, model.row_lb, model.row_ub, lbNew, ubNew);
%ResMaxBiomass=tomRun('cplex', ProbMaxBiomass, 0);
ResMaxBiomass=runLP(model.organismObjective, model.S, model.row_lb, model.row_ub, lbNew, ubNew);


ubNew(model.organismObjectiveInd)=ResMaxBiomass.obj;
lbNew(model.organismObjectiveInd)=ResMaxBiomass.obj;

if (ResMaxBiomass.stat == 1)
    %res.chemical = model.chemicalInd;
    res.knockedInds = knockedInds;
    %res.maxW=maxW;
    res.objectiveVal=ResMaxBiomass.obj;

    %maximum secretion and matching fluxes
    %Prob_maxChemical=lpAssign(-model.C_chemical, model.S, model.row_lb, model.row_ub, lbNew, ubNew);
    %Result_maxChemical=tomRun('cplex', Prob_maxChemical, 0);
    Result_maxChemical=runLP(model.C_chemical, model.S, model.row_lb, model.row_ub, lbNew, ubNew);

    %res.exitFlagOfMaxChemical=Result_maxChemical.ExitFlag;
    res.maxChemical=Result_maxChemical.full(model.chemicalInd);
    %res.x_kOfMaxChemical=Result_maxChemical.x_k;

    %minimum secretion and fluxes
    %Prob_minChemical=lpAssign(model.C_chemical, model.S, model.row_lb, model.row_ub, lbNew, ubNew);
    %Result_minChemical=tomRun('cplex', Prob_minChemical, 0);
    %res.exitFlagOfGuerenteedChemical=Result_minChemical.ExitFlag;
    Result_minChemical=runLP(-model.C_chemical, model.S, model.row_lb, model.row_ub, lbNew, ubNew);
    res.GuaranteedChemical=Result_minChemical.full(model.chemicalInd);
    %res.x_kOfMinChemical=Result_minChemical.x_k;
end
