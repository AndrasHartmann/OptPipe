
% Read the model, insert heterologous reactions and constraints
model = readCbModel([cd filesep 'Glutamicum_models' filesep 'ReliableModel_ElisabethZelle' filesep 'iEZ475.xml']);
[cobra_model tmp]=Get_heterologous(model, 'pelargonidin');
sink = 'sink_naringenin';

% parameter settings
BM_rxn='biomass_a';
target_rxn= sink;

cobra_model=Insert_constraints(cobra_model);
cobra_model=changeRxnBounds(cobra_model,'sink_coumarate',-0.5,'l');

productionEnvelope(cobra_model, {},'k--',sink,BM_rxn)
hold on

productionEnvelope(cobra_model, {'sdhCAB','ddh'},'r--',sink,BM_rxn)
productionEnvelope(cobra_model, {'sdhCAB'},'b--',sink,BM_rxn)
productionEnvelope(cobra_model, {'ddh'},'g--',sink,BM_rxn)

legend('wild-type', 'sdhCAB-ddh', 'sdhCAB', 'ddh')
