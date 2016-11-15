function writeRankedReactions(cobra_model, ranked_reactions, filename, ko)

    fileID = fopen(filename, 'w');
    
    if ko
        % Knockout:
        fprintf(fileID, '%s\t%s\t%s\t%s\t%s\t%s\t%s\n', 'Reaction id', 'Reaction name', 'Reaction definition', 'Gene association', 'High flux loss?', 'KO ranking', 'Subsystem');
    else
        % Overexpression:
        fprintf(fileID, '%s\t%s\t%s\t%s\t%s\t%s\n', 'Reaction id', 'Reaction name', 'Reaction definition', 'Gene association', 'Overexpression ranking', 'Subsystem');
    end
    
    for row = 1:size(ranked_reactions, 1)
        reaction_index = ranked_reactions(row, 1);
        fprintf(fileID, '%s\t%s\t%s\t%s\t', char(cobra_model.rxns(reaction_index)), char(cobra_model.rxnNames(reaction_index)), char(printRxnFormula(cobra_model, cobra_model.rxns(reaction_index), false)), char(cobra_model.grRules(reaction_index)));
        
        if ko
            fprintf(fileID, '%d\t%f\t', ranked_reactions(row, 2), ranked_reactions(row, 3));
        else
            fprintf(fileID, '%f\t', ranked_reactions(row, 2));
        end
        
        fprintf(fileID, '%s\n', char(cobra_model.subSystems(reaction_index)));
    end
    
    fclose(fileID);
end