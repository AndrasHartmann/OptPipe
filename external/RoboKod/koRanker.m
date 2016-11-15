function [ranked_reactions, range, OTDB_min, OTDB_max] ...
    = koRanker(cobra_model, biomass_reaction_id, target_reaction_id, important_reactions)

cobra_model = changeObjective(cobra_model, biomass_reaction_id);

% Ignore exchange reactions:
exchange_reaction_ids = find(findExcRxns(cobra_model,1));
non_exchange_reaction_ids = setdiff(1:length(cobra_model.rxns), exchange_reaction_ids);

% Optimise target while making biomass at different percentages (from 0 to 100%) of maximum...
range = 0:5:100;
[OTDB_min, OTDB_max] = multiFluxVariability(cobra_model, non_exchange_reaction_ids, biomass_reaction_id, target_reaction_id, range);

ranked_reactions = [];

number_fva_points = size(OTDB_min, 2);
percentile_5 = ceil(number_fva_points * 5/100);
percentile_10 = ceil(number_fva_points * 10/100);
percentile_25 = ceil(number_fva_points * 25/100);
percentile_75 = ceil(number_fva_points * 75/100);

for k = 1:length(non_exchange_reaction_ids)
    % Exclude reactions vital for target production (i.e. those that have
    % a non-zero max or min flux at high target production):
    target_zero = [OTDB_max(k, percentile_5:percentile_10); OTDB_min(k, percentile_5:percentile_10)];

    if any(abs(target_zero) < 1e-6)
        OTDB_min_percentile_25 = OTDB_min(k, 1:percentile_25);
        OTDB_max_percentile_25 = OTDB_max(k, 1:percentile_25);
        OTDB_min_percentile_75 = OTDB_min(k, percentile_75:end);
        OTDB_max_percentile_75 = OTDB_max(k, percentile_75:end);
        
        if sum(OTDB_min_percentile_25) < 0 && sum(OTDB_max_percentile_25) <= 0
            target_ex_ub_diff = sum(abs(OTDB_min_percentile_75) - abs(OTDB_min_percentile_25));
            area_ub = sum(abs(OTDB_min_percentile_75) - abs(OTDB_max_percentile_75));
        elseif sum(OTDB_min_percentile_25) >= 0 && sum(OTDB_max_percentile_25) > 0
            target_ex_ub_diff = sum(OTDB_max_percentile_75 - OTDB_max_percentile_25);
            area_ub = sum(OTDB_max_percentile_75 - OTDB_min_percentile_75);
        else
            continue
        end
        
        computeTotal = target_ex_ub_diff / area_ub;
        
        if computeTotal > 0 && ~isinf(computeTotal)
        	reaction_index = non_exchange_reaction_ids(k);
        	ranked_reactions = [ranked_reactions; reaction_index, any(important_reactions == reaction_index), computeTotal];
        end
    end
end

% Sort:
if ~isempty(ranked_reactions)
    ranked_reactions = sortrows(ranked_reactions, [-3 -2]);
    
    % Only those FVA values for ranked_reactions need to be returned.
    % Map indexes of ranked_reactions to non_exchange_reaction_ids:
    fva_reaction_indices = arrayfun(@(x) find(non_exchange_reaction_ids == x, 1), ranked_reactions(:,1));
    OTDB_min = OTDB_min(fva_reaction_indices, :);
    OTDB_max = OTDB_max(fva_reaction_indices, :);
else
    OTDB_min = [];
    OTDB_max = [];
end