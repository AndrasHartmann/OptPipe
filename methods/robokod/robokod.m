function [results] = robokod(cobra_model, biomass_reaction_id, target_reaction_id, max_knockouts,selectedRxnIDs)
    %Init
    results={};
    equivalent_KOs='';
    
    %% Iterative knockout analysis:
    for i = 1:max_knockouts
        
        disp(strcat('Initiating KO iteration number ',num2str(i)));

        % MCT:
        important_reactions = metaboliteCentralityTest(cobra_model, biomass_reaction_id, target_reaction_id);

        % FVAp:
        [ranked_reactions, range, OTDB_fva_min, OTDB_biomass_fva_max] ...
            = koRanker(cobra_model, biomass_reaction_id, target_reaction_id, important_reactions,selectedRxnIDs);

        % Exit iterative loop if no candidate knock-out reactions returned:
        if isempty(ranked_reactions)
           break; 
        end
        
            
        %% Check if there are equivalent KOs: save the equivalencies
        %find the 1st adjacent KOs with different rankings
        big_diffs=find(abs(diff(ranked_reactions(:,3)))>0.1);

        %In case the highest ranked KO has equivalent KOs,
        %find the equivalencies and save them
        if big_diffs(1)~=1
            equivs=cobra_model.rxns(ranked_reactions((1:big_diffs(1)),1));
            %equivalent_KOs{end+1,1}=strjoin(equivs',',');
            %Just to guarantee compatibility with R2012a
            %TODO: check
            joinedstr = sprintf('%s,' ,equivs{:});
            joinedstr(end) = [];
            equivalent_KOs{end+1,1} = joinedstr;
        end
        
        %% Save results and prepare KO
        results=horzcat(results,cobra_model.rxns(ranked_reactions(1,1)));
 
        % Knockout the highest ranked reaction:
        cobra_model = changeRxnBounds(cobra_model, cobra_model.rxns(ranked_reactions(1,1)), 0, 'b');
        
    end
    
   
    %% Organize the results: create the possible combinations with equivalent KOs
    combination=results(1,:);

    for i=1:length(equivalent_KOs)
        pair=equivalent_KOs{i,1};
        %KOs=strsplit(pair,',');    %the 1st KO is the highest ranked
        %Just to guarantee compatibility with R2012a
        %TODO: check
        KOs = regexp (pair, ',', 'split');
        for j=2:length(KOs)
            new_comb=combination;
            new_comb( strmatch(KOs{1,1},combination) ) = KOs(j);
            results(end+1,:)=new_comb;
        end
    end



end
