%Script to generate the Envelopes
changeCobraSolver('gurobi5');

% Read the model, insert heterologous reactions and constraints
model = readCbModel([cd filesep 'Glutamicum_models' filesep 'ReliableModel_ElisabethZelle' filesep 'iEZ475.xml']);
[cobra_model tmp]=Get_heterologous(model, 'pelargonidin');
sink = 'sink_naringenin';

% parameter settings
BM_rxn='biomass_a';
target_rxn= sink;

optpipe_model=Insert_constraints(cobra_model);
optpipe_model=changeRxnBounds(optpipe_model,'sink_coumarate',-0.5,'l');

figure
productionEnvelope(optpipe_model, {},'k--',sink,BM_rxn)
hold on

productionEnvelope(optpipe_model, {'sdhCAB','ddh'},'r--',sink,BM_rxn)
productionEnvelope(optpipe_model, {'sdhCAB'},'b--',sink,BM_rxn)
productionEnvelope(optpipe_model, {'ddh'},'g--',sink,BM_rxn)

title ('Naringenin with coumarate and constraints');
legend('wild-type', 'sdhCAB-ddh', 'sdhCAB', 'ddh')


%% Malonyl-coA model
mal_model=Insert_constraints(model);
mal_model=addReaction(mal_model,'EX_malcoa',{'malcoa[c]','malcoa[e]'},[-1 1]);
mal_model.rxnNames{end,1}='EX malcoa';
mal_model=changeRxnBounds(mal_model,'EX_malcoa',0,'l');
mal_model=addReaction(mal_model,'sink_malcoa',{'malcoa[e]'},[-1]);
mal_model.rxnNames{end,1}='sink malcoa'


sink='sink_malcoa';

hold off
figure
productionEnvelope(mal_model, {},'k--',sink,BM_rxn)
hold on
productionEnvelope(mal_model, {'sdhCAB','ddh'},'r--',sink,BM_rxn)
productionEnvelope(mal_model, {'sdhCAB'},'b--',sink,BM_rxn)
productionEnvelope(mal_model, {'ddh'},'g--',sink,BM_rxn)
productionEnvelope(mal_model, {'mdh','glgA','aspA'},'m--',sink,BM_rxn)

title ('Malonyl coA');
legend('wild-type', 'sdhCAB-ddh', 'sdhCAB', 'ddh','mdh-glgA-aspA')





% 
% %%No coumrate (noCoumarate_model)
% 
% noCoumarate_model=Insert_constraints(cobra_model);
% 
% 
% 
% hold off
% subplot(3,1,2);
% productionEnvelope(noCoumarate_model, {},'k--',sink,BM_rxn)
% hold on
% 
% productionEnvelope(noCoumarate_model, {'sdhCAB','ddh'},'r--',sink,BM_rxn)
% productionEnvelope(noCoumarate_model, {'sdhCAB'},'b--',sink,BM_rxn)
% productionEnvelope(noCoumarate_model, {'ddh'},'g--',sink,BM_rxn)
% 
% title ('No coumarate');
% legend('wild-type', 'sdhCAB-ddh', 'sdhCAB', 'ddh')
% 
% %not this
% 
% %% %Check rigid contraints and with coumarate (orig_model)
% 
% orig_model = Insert_constraints( cobra_model,true)
% orig_model=changeRxnBounds(orig_model,'sink_coumarate',-0.5,'l');
% 
% hold off
% subplot(3,1,3);
% productionEnvelope(orig_model, {},'k--',sink,BM_rxn)
% hold on
% 
% productionEnvelope(orig_model, {'sdhCAB','ddh'},'r--',sink,BM_rxn)
% productionEnvelope(orig_model, {'sdhCAB'},'b--',sink,BM_rxn)
% productionEnvelope(orig_model, {'ddh'},'g--',sink,BM_rxn)
% 
% title ('With coumarate and rigid constraints');
% legend('wild-type', 'sdhCAB-ddh', 'sdhCAB', 'ddh')
% %not this
