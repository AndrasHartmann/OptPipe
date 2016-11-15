function [biomass target_min target_max distance] = validateDeletion(model, delRxns, targetRxn,FVA)

%addpath('../enumerate');

% Setting the objective for target
p = zeros(size(model.c));
targetRxnInd = find(ismember(model.rxns,  targetRxn));
p(targetRxnInd) = -1;
delmodel = delReaction(model, find(ismember(model.rxns, delRxns)));

tmp = robustKnockSolution(delmodel,p);

biomass = tmp(1);
target_max = tmp(2);

p(targetRxnInd) = 1;

tmp = robustKnockSolution(delmodel,p);

target_min = tmp(2);


%Compute distance between Wt and mutant
FBA=optimizeCbModel(delmodel);
distance=0;
for i=1:length(FVA.min) 
    dif=0;
    flux=FBA.x(i);
    if flux < FVA.min(i)
        if sign(flux)==sign(FVA.min(i))
            dif=abs(flux)-abs(FVA.min(i));
        else
            dif=FVA.min(i)-flux;
        end
    end
    if flux > FVA.max(i)
        if sign(flux)==sign(FVA.max(i))
            dif=abs(flux)-abs(FVA.max(i));
        else
            dif=flux-FVA.max(i);
        end
    end
    %only take into account when the dif is over 10% of minimum
    tol=0.1*abs(FVA.min(i));
    if dif>tol
        distance=distance+dif;
    end
end
end


%TODO: if we want to plot figure
%{
set(0,'DefaultAxesFontSize', 15)
set(0,'DefaultAxesFontWeight', 'bold')
set(0,'DefaultTextFontSize', 15)
set(0,'DefaultTextFontWeight', 'bold')

figure
hold on;
[a, b] = productionEnvelope(model,{}, 'b--',model.rxns(end),model.rxns(301),false,20)
[a, b] = productionEnvelope(model,{'sdhCAB',     'ddh'     }, 'r--',model.rxns(end),model.rxns(301),false,20)
[a, b] = productionEnvelope(model,{'mdh',    'glgC',    'aspA' }, 'k',model.rxns(end),model.rxns(301),false,20)
[a, b] = productionEnvelope(model,{'mdh',    'glgA',    'aspA' }, 'g-.',model.rxns(end),model.rxns(301),false,20)
[a, b] = productionEnvelope(model,{'mdh',    'aspA',    'NH3_d'}, 'm--',model.rxns(end),model.rxns(301),false,20)

l = legend( 'Wild type', '$\Delta$sdhCAB$\Delta$ddh', '$\Delta$mdh$\Delta$glgA$\Delta$aspA', '$\Delta$mdh$\Delta$glgC$\Delta$aspA', '$\Delta$mdh$\Delta$NH3\_d$\Delta$aspA') 
set(l,'Interpreter','latex');

xlabel('growth rate (h^-1)')
ylabel('McoA (mmol/gDW/h)')
%}
