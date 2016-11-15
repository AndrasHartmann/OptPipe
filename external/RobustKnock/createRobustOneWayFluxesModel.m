function [model,yInd, notYInd,  m,n, coupled, coupledYInd, coupledNotYInd]= ...
createRobustOneWayFluxesModel(originalModel, chemicalInd, coupledFlag, knockablerxns)
%function [model,yInd, notYInd,  m,n, coupled, coupledYInd, coupledNotYInd]= createRobustOneWayFluxesModel(originalModel, chemicalInd, coupledFlag, knockablerxns)
% INPUTS
% originalModel
% chemicalInd
% coupledFlag - 1 for a reversible model; 0 for an
%               irreversible model.
%
% OUTPUTS
% originalModel
% yInd - Can be knocked out
% notYInd - Cannot be knocked out
% m - size(model.S,1)
% n - size(model.S,2)
% coupled - Coupled to GPR
% coupledYInd - Coupled to GPR, can be knocked out
% coupledNotYInd - Coupled to GPR, cannot be knocked out

%% Initialization

tmpModel=originalModel;
[met_num rxn_num] = size(tmpModel.S);
yInd=[];
notYInd=[];
coupled=[];
coupledYInd=[];
coupledNotYInd=[];
coupledInd=[];
%vEqual=[];
%seperateFlux=[];


%%Construct one directional model
%create variables>=0
%matrices to change: S,c,lb,ub,rxns,rxnGeneMat, rxnNames,  in tmpModel
tmp_rxn_num=rxn_num;
if (coupledFlag)
    for i=1:rxn_num
        lbound=tmpModel.lb(i);
        ubound=tmpModel.ub(i);
        %if lower bound smaller than zero
        if(lbound < 0)
            sLeft=tmpModel.S(:,1:i-1);
            sRight=tmpModel.S(:,i+1:tmp_rxn_num);
            sVar=tmpModel.S(:,i);

            rxnGeneMatLeft=tmpModel.rxnGeneMat(1:i-1,:);
            rxnGeneMatRight=tmpModel.rxnGeneMat(i+1:tmp_rxn_num,:);
            rxnGeneMatVar=tmpModel.rxnGeneMat(i,:);

            cLeft=tmpModel.c(1:i-1);
            cRight=tmpModel.c(i+1:tmp_rxn_num);
            cVar=tmpModel.c(i);
            
            organismObjectiveLeft=tmpModel.organismObjective(1:i-1);
            organismObjectiveRight=tmpModel.organismObjective(i+1:tmp_rxn_num);
            organismObjectiveVar=tmpModel.organismObjective(i);
            
            C_chemicalLeft=tmpModel.C_chemical(1:i-1);
            C_chemicalRight=tmpModel.C_chemical(i+1:tmp_rxn_num);
            C_chemicalVar=tmpModel.C_chemical(i);
       
            lbLeft=tmpModel.lb(1:i-1);
            lbRight=tmpModel.lb(i+1:tmp_rxn_num);

            ubLeft=tmpModel.ub(1:i-1);
            ubRight=tmpModel.ub(i+1:tmp_rxn_num);

            rxnsLeft=tmpModel.rxns(1:i-1);
            rxnsRight=tmpModel.rxns(i+1:tmp_rxn_num);
            rxnsVar=tmpModel.rxns(i);

            rxnNamesLeft=tmpModel.rxnNames(1:i-1);
            rxnNamesRight=tmpModel.rxnNames(i+1:tmp_rxn_num);
            rxnNamesVar=tmpModel.rxnNames(i);

            if(ubound>0)    %need to create 2 reactions (split the reaction in two)
                if (i == chemicalInd)
                    tmpModel.lb(i)=0;  % only secretion from this chemical!
                else
                    s1=sVar;
                    s2=-sVar;
                    c1=cVar;
                    c2=-cVar;
                    lb1=0;
                    lb2=0;
                    ub1=ubound;
                    ub2=-lbound;
                    organismObjective1=organismObjectiveVar;
                    organismObjective2=-organismObjectiveVar;
                    C_chemical1=C_chemicalVar;
                    C_chemical2=-C_chemicalVar;                    

                    Snew=[sLeft,s1,sRight,s2];
                    cnew=[cLeft;c1;cRight;c2];
                    lbnew=[lbLeft;lb1;lbRight;lb2];
                    ubnew=[ubLeft;ub1;ubRight;ub2];
                    rxnsNew=[rxnsLeft;rxnsVar;rxnsRight;strcat(rxnsVar,'2')];
                    rxnNamesNew=[rxnNamesLeft;rxnNamesVar;rxnNamesRight;strcat(rxnNamesVar,'2')];
                    rxnGeneMatNew=[rxnGeneMatLeft;rxnGeneMatVar;rxnGeneMatRight;rxnGeneMatVar];
                    organismObjectiveNew=[organismObjectiveLeft;organismObjective1;organismObjectiveRight;organismObjective2];
                    C_chemicalNew=[C_chemicalLeft;C_chemical1;C_chemicalRight;C_chemical2];

                    tmpModel.S=Snew;
                    tmpModel.c=cnew;
                    tmpModel.lb=lbnew;
                    tmpModel.ub=ubnew;
                    tmpModel.rxns=rxnsNew;
                    tmpModel.rxnNames=rxnNamesNew;
                    tmpModel.rxnGeneMat=rxnGeneMatNew;
                    tmpModel.organismObjective=organismObjectiveNew;
                    tmpModel.C_chemical=C_chemicalNew;

                    tmp_rxn_num=tmp_rxn_num+1;
                    %couple the original and the new reaction
                    %e.g. they can be cknocked out together
                    coupled(end+1, 1:2)=[i,tmp_rxn_num];
                end
            end
            if(ubound<=0)    %need to change directions
                Snew=[sLeft,-sVar,sRight];
                cnew=[cLeft;-cVar;cRight];
                lbnew=[lbLeft;-ubound;lbRight];
                ubnew=[ubLeft;-lbound;ubRight];
                organismObjectiveNew=[organismObjectiveLeft;-organismObjectiveVar;organismObjectiveRight];
                C_chemicalNew=[C_chemicalLeft;C_chemicalVar;C_chemicalRight];
                
                tmpModel.S=Snew;
                tmpModel.c=cnew;
                tmpModel.lb=lbnew;
                tmpModel.ub=ubnew;
                tmpModel.organismObjective=organismObjectiveNew;
                tmpModel.C_chemical=C_chemicalNew;
            end
        end
    end
end
% construct constraints
model = tmpModel;
[n m] = size(model.S);

for i=1:rxn_num
    if (find(notYInd == i) )
        %go to the next index
        continue;
    else
        [xCoupled,y]=find(coupled==i);
        if(xCoupled)
            coupledInd=coupled(xCoupled, y+1);
        end
        if find(knockablerxns == i)
            yInd(end+1)=i;
            if(xCoupled)
                coupledYInd(end+1)=coupledInd;
            end
            continue;
        else
            notYInd(end+1)=i;
            if(xCoupled)
                coupledNotYInd(end+1)=coupledInd;
            end
        end
    end
end

notYInd=sort(notYInd)';
yInd=sort(yInd)';
coupledNotYInd=sort(coupledNotYInd)';
coupledYInd=sort(coupledYInd)';

%Sorting
%yInd = sort(yInd)'; yCoupledInd = sort(yCoupledInd)';
%qInd = sort(qInd)'; qCoupledInd = sort(qCoupledInd)';
%sInd = sort(sInd)'; sCoupledInd = sort(sCoupledInd)';
