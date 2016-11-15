function [importantReactions] = metaboliteCentralityTest(cobra_model, biomass_reaction_id, target_reaction_id)

% need the FBA result of the model when generating biomass, so that the
% flux through each reaction can be identified.
cobra_model = changeObjective(cobra_model, biomass_reaction_id);
biomass_objective_solution = optimizeCbModel(cobra_model, 'max');
[min,max] = fluxVariability(cobra_model,95,'max');
zeroFluxRxns = (biomass_objective_solution.x == 0);

%---------------------------------------
%Identify active target reactions
%---------------------------------------
cobra_model = changeObjective(cobra_model, target_reaction_id);
target_solution = optimizeCbModel(cobra_model, 'max', 'one');
target_reactions = abs(target_solution.x) > 1e-6;

%---------------------------------------
% cycle through the reactions. Identify
% metabolites made in each reaction
% See whether target reactions use the
% metabolites
%---------------------------------------

N = full(cobra_model.S);
newN = (N(:,target_reactions));
t = sum(abs(newN'));
all_species_involved = find(t ~= 0);

% Remove cofactors:
bounds = [{'adp[c]'}, {'amp[c]'}, {'atp[c]'}, {'co2[c]'}, {'co2[e]'},...
    {'coa[c]'}, {'dhap[c]'}, {'e4p[c]'}, {'etoh[e]'}, {'h2o[c]'},...
    {'h2o[e]'}, {'h[c]'}, {'h[e]'}, {'nad[c]'}, {'nadh[c]'},...
    {'nadp[c]'}, {'nadph[c]'}, {'nh4[e]'}, {'o2[c]'}, {'o2[e]'},...
    {'pi[c]'}, {'pi[e]'}];

bounds_metabolite_indices = find(ismember(cobra_model.mets, bounds));
species_involved = setdiff(all_species_involved, bounds_metabolite_indices);

sums = [{},{},{}];

for k = 1:length(species_involved)
    % count up flux loss
    
    count = 0;
    fluxCount = 0;
    
    %need to make sure you are getting all reactions.
    reactions_involved = find(N(species_involved(k),:)~=0);
    
    for kk = 1:length(reactions_involved)

        target_reaction_flux = target_solution.x(reactions_involved(kk));
        biomass_reaction_flux = biomass_objective_solution.x(reactions_involved(kk));
        variability = zeroFluxRxns(reactions_involved(kk));
        maxFVA = max(reactions_involved(kk));

        if(abs(target_reaction_flux) > 1e-6) && (abs(biomass_reaction_flux) > 1e-6)
            if(biomass_reaction_flux*target_reaction_flux > 0)
                %---------------------------------------
                % Both target and biomass carry flux
                % in the same direction so no flux lost.
                %---------------------------------------
            else
                %---------------------------------------
                % Both biomass and target carry flux
                % but in opposite directions.
                %---------------------------------------
                continue
            end

        %---------------------------------------
        % target doesn't use reaction, biomass does in
        % a reverse direction and reaction not needed.
        %---------------------------------------
        elseif (abs(target_reaction_flux) < 1e-6) && biomass_reaction_flux < 0 && (abs(variability) == 1)
            if (N(species_involved(k),reactions_involved(kk)) < 0) &&  (abs(biomass_reaction_flux) >= 1)
                count = count+1;
                fluxCount = fluxCount+(count*abs(maxFVA));
            elseif (N(species_involved(k),reactions_involved(kk)) < 0) &&  (abs(biomass_reaction_flux) < 1) &&  (abs(biomass_reaction_flux) > 1e-4)
                count = count+1;
                fluxCount = fluxCount+(count*abs(maxFVA));
            elseif (N(species_involved(k),reactions_involved(kk)) > 0) && (abs(biomass_reaction_flux) >= 1)
                count = count-1;
                fluxCount = fluxCount+(count*abs(maxFVA));
            elseif (N(species_involved(k),reactions_involved(kk)) > 0) &&  (abs(biomass_reaction_flux) > 1e-4) &&  (abs(biomass_reaction_flux) < 1)
                count = count-1;
                fluxCount = fluxCount+(count*abs(maxFVA));
            end
        %---------------------------------------
        % target doesn't use reaction, biomass does in
        % a forwards direction, and the variability
        % means reaction is MAYBE not needed.
        %---------------------------------------
        elseif (abs(target_reaction_flux) < 1e-6) && (abs(variability) == 1)
            if (N(species_involved(k),reactions_involved(kk)) < 0) && (biomass_reaction_flux > 1e-4) && (biomass_reaction_flux < 1)
                count = count-1;
                fluxCount = fluxCount+(count*abs(maxFVA));
            elseif (N(species_involved(k),reactions_involved(kk)) < 0) && (biomass_reaction_flux >= 1)
                count = count-1;
                fluxCount = fluxCount+(count*abs(maxFVA));
            elseif (N(species_involved(k),reactions_involved(kk)) > 0) && (biomass_reaction_flux >= 1)
                count = count+1;
                fluxCount = fluxCount+(count*abs(maxFVA));
            elseif (N(species_involved(k),reactions_involved(kk)) > 0) && (biomass_reaction_flux > 1e-4) && (biomass_reaction_flux < 1)
                count = count+1;
                fluxCount = fluxCount+(count*abs(maxFVA));
            end
        %---------------------------------------
        % target doesn't use reaction, biomass does in
        % a forwards direction, and the variability
        % means reaction is needed.
        %---------------------------------------
        elseif (abs(target_reaction_flux) < 1e-6) && (biomass_reaction_flux > 1e-6) && (abs(variability) == 0)
            if (N(species_involved(k),reactions_involved(kk)) < 0) && (biomass_reaction_flux >= 1)
                count = count-1;
                fluxCount = fluxCount+(count*abs(maxFVA));
            elseif (N(species_involved(k),reactions_involved(kk)) < 0) && (biomass_reaction_flux > 1e-4) && (biomass_reaction_flux < 1)
                count = count-1;
                fluxCount = fluxCount+(count*abs(maxFVA));
            elseif (N(species_involved(k),reactions_involved(kk)) > 0) && (biomass_reaction_flux >= 1)
                count = count+1;
                fluxCount = fluxCount+(count*abs(maxFVA));
            elseif (N(species_involved(k),reactions_involved(kk)) > 0) && (biomass_reaction_flux < 1) && (biomass_reaction_flux > 1e-4)
                count = count+1;
                fluxCount = fluxCount+(count*abs(maxFVA));
            end
        %---------------------------------------
        % target doesn't use reaction, biomass does in
        % a reverse direction and variability is there.
        %---------------------------------------
        elseif (abs(target_reaction_flux) < 1e-6) && (abs(biomass_reaction_flux) > 1e-6) && (biomass_reaction_flux < 0) && (abs(variability) == 0)
            if (N(species_involved(k),reactions_involved(kk)) < 0) && (abs(biomass_reaction_flux) >= 1)
                count = count+1;
                fluxCount = fluxCount+(count*abs(maxFVA));
            elseif (N(species_involved(k),reactions_involved(kk)) < 0) && (abs(biomass_reaction_flux) < 1) && (abs(biomass_reaction_flux) > 1e-4)
                count = count+1;
                fluxCount = fluxCount+(count*abs(maxFVA));
            elseif (N(species_involved(k),reactions_involved(kk)) > 0) && (abs(biomass_reaction_flux) >= 1)
                count = count-1;
                fluxCount = fluxCount+(count*abs(maxFVA));
            elseif (N(species_involved(k),reactions_involved(kk)) > 0) && (abs(biomass_reaction_flux) < 1) && (abs(biomass_reaction_flux) > 1e-4)
                count = count-1;
                fluxCount = fluxCount+(count*abs(maxFVA));
            end
        end
    end
    
    sums = [sums; {cobra_model.mets(species_involved(k))},{fluxCount}];
    sums = sortrows(sums,2);
end

%---------------------------------------
% output
%---------------------------------------
% select the highest flux loss metabolite and identify all reactions
% around it.
metabolite_index = strcmp(cobra_model.mets, sums{1,1});
[~,importantReactions] = find(cobra_model.S(metabolite_index, :));