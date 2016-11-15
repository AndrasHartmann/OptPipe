function [selectedRxns] = getSelectedRxns(model,biomassRxn, onlyTrans)
%function [selectedRxns] = getSelectedRxns(model,biomassRxn, onlyTrans)
%
%This function is for pre-processing. It determines a list of "selected" reactions that should not be considered for Knock-Out.
% Irrelevant reactions that are not contained in selected reactions are the essential, transport and blocked blocked reactions. 

%% Initialization
%Just in case objective would be different
cobra_model=changeObjective(model,biomassRxn);

%% If only the tranport reactions are needed
if exist('onlyTrans')&&onlyTrans
    selectedRxns = cobra_model.rxns(~ismember(cobra_model.rxns, findTransRxns(model,true)));
    return;
end



%%Find biomass rxn
biomassrxn = (ismember(model.rxns,biomassRxn));
% lower flux limit for reaction to be considered lethal 
fluxTol = 1e-9; 

%% Lethal Rxns
[~,~,~,~,~,Flux_FBA] = singleRxnDeletion(model,'FBA');
lethalRxns = Flux_FBA( biomassrxn,: )< fluxTol; 

%% Transport and synthetic reactions
transRxns = ismember(model.rxns, findTransRxns(cobra_model,true));

%% Transport and synthetic reactions
excRxns = findExcRxns(cobra_model,true);
    
%% Blocked reactions
blockedRxns = ismember(model.rxns, findBlockedReaction(model))';

%%join the results
irrelevant = excRxns(:)|blockedRxns(:)|transRxns(:)|lethalRxns(:)|biomassrxn(:);
selectedRxns = model.rxns(~irrelevant);

