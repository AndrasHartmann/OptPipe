function robokod(model_filename, biomass_reaction_id, target_reaction_id, max_knockouts, results_directory)
    % Make results directory:
    if ~exist(fullfile(cd, results_directory), 'file')
        mkdir(results_directory);
    end

    % Read model:
    cobra_model = readCbModel(model_filename);

    %% Iterative knockout analysis:
    for i = 1:max_knockouts

        % MCT:
        important_reactions = metaboliteCentralityTest(cobra_model, biomass_reaction_id, target_reaction_id);

        % FVAp:
        [ranked_reactions, range, OTDB_fva_min, OTDB_biomass_fva_max] ...
            = koRanker(cobra_model, biomass_reaction_id, target_reaction_id, important_reactions);

        % Exit iterative loop if no candidate knock-out reactions returned:
        if isempty(ranked_reactions)
           break; 
        end

        % Write results to Excel / plot FVAp distributions:
        writeRankedReactions(cobra_model, ranked_reactions, strcat(results_directory, '/', 'ko_ranked_reactions_', num2str(i), '.xls'), 1);
        fva_plot_directory_name = strcat('fva_plots_', num2str(i));
        fva_plot_directory_path = strcat(results_directory, '/', fva_plot_directory_name);

        if ~exist(fullfile(cd, fva_plot_directory_path), 'file')
            mkdir(results_directory, fva_plot_directory_name);
        end

        % Optimise biomass while making target at different percentages (from 0 to % 99%) of maximum...
        % (This is purely for plotting).
        [OBDT_fva_min, OBDT_fva_max] = multiFluxVariability(cobra_model, ranked_reactions(:,1), target_reaction_id, biomass_reaction_id, range);
        OBDT_fva_min = fliplr(OBDT_fva_min);
        OBDT_fva_max = fliplr(OBDT_fva_max);
        fva_plot(ranked_reactions(:,1), range, OTDB_biomass_fva_max, OTDB_fva_min, OBDT_fva_max, OBDT_fva_min, cobra_model, fva_plot_directory_path);

        % Knockout the highest ranked reaction:
        cobra_model = changeRxnBounds(cobra_model, cobra_model.rxns(ranked_reactions(1,1)), 0, 'b');

        % Output to console:
        fprintf( 'Knockout %d: %s (%s)\n', i, char(cobra_model.rxnNames(ranked_reactions(1,1))), char(cobra_model.rxns(ranked_reactions(1,1))));
    end

    % Write updated SBML model including knockouts:
    if max_knockouts > 0
        model_output_filepath = strcat(results_directory, '/', strrep(model_filename, '.xml', '_ko.xml'));
        writeCbToSBML(cobra_model, model_output_filepath);
    end


    %% Calculate fluxes on updated model
    model = changeObjective(cobra_model, biomass_reaction_id);

    FBAsolution = optimizeCbModel(model, 'max', 'one');

    if FBAsolution.f == 0
        disp('No solution found')
    else
        result_output_filepath = strcat(results_directory, '/', strrep(model_filename, '.xml', '_ko_biomass_fba_results.xls'));
        writeResult(model, FBAsolution, result_output_filepath);
    end


    %% Overexpression analysis
    [ranked_reactions, range, OTDB_fva_min, OTDB_biomass_fva_max, OBDT_fva_min, OBDT_fva_max] ...
        = overexpressRanker(cobra_model, biomass_reaction_id, target_reaction_id);

    % Write results to Excel / plot FVAp distributions:
    writeRankedReactions(cobra_model, ranked_reactions, strcat(results_directory, '/', 'overexpression_ranked_reactions.xls'), 0);

    fva_plot_directory_name = 'fva_plots_overexpression';
    fva_plot_directory_path = strcat(results_directory, '/', fva_plot_directory_name);

    if ~exist(fullfile(cd, fva_plot_directory_path), 'file')
        mkdir(results_directory, fva_plot_directory_name);
    end

    fva_plot(ranked_reactions(:,1), range, OTDB_biomass_fva_max, OTDB_fva_min, OBDT_fva_max, OBDT_fva_min, cobra_model, fva_plot_directory_path);

    % Output to console:
    fprintf( 'Overexpression: %d candidate reactions\n', length(ranked_reactions));
end