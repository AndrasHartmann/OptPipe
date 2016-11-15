function [delRxn,fluxSolution] = testRxnDeletion(model,evaluate,numDel,rxnList)
%function [delRxn,fluxSolution] = testRxnDeletion(model,evaluate,numDel,rxnList)

if (nargin < 4)
    error('testRxnDeletion needs at least 3 arguments')
end

if (nargin < 4)
    rxnList = 1:(size(model.rxns));
else
    if (isempty(rxnList))
        rxnList = 1:(size(model.rxns));
    end
end

tol = 10^-9;

nDelRxns = length(rxnList);
    if numDel == 1
        delRxn = zeros(nDelRxns,1);
        fluxSolution = zeros(nDelRxns,2);

        parfor i = 1:nDelRxns
            changeCobraSolver('gurobi6');
            delRxn(i) = rxnList(i);
            fluxSolution(i,:) = evaluate(delReaction(model, rxnList(i)));
        end;
        return;
    else
        delRxn = [];
        fluxSolution = [];
        for i = 1:nDelRxns-numDel+1
            %Uncomment this if you would like to have some status
            %if numDel == 3
            %    i
            %end
            [delRxn_,fluxSolution_] = testRxnDeletion(delReaction(model, rxnList(i)), evaluate, numDel-1, rxnList(i+1:end));
            
            delRxn_=[ones(size(delRxn_,1),1)*rxnList(i) delRxn_];
            %Only those solutions where biomass larger than tolerance
            delRxn = [delRxn; delRxn_( find(fluxSolution_(:,2)>tol), :) ];
            fluxSolution = [fluxSolution; fluxSolution_( find(fluxSolution_(:,2)>tol), :)];
        end;
    end;

end

