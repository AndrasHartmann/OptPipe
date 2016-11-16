function P=rankoptpipe(M)

%RANKOPTPIPE calculates rank product of meta-optimization studies
% function P=rankoptpipe(M)
% INPUT
%   M - matrix (n x k), each row is deletion, column is
%   (experiments/methods), higher values are better, e.g. highest value in
%   column gets rank=1). Pre-process column with symmetric values if lowest
%   is better.
% OUTPUT
%   P - struct with p-value and q-value information
%       .PG - gamma aproximation
%       .PE - exact calculation
%       .QG - q-value of PG
%       .QE - q-value of QE
% Calls external function -  piltzcount.m
% See also: OptPipe toolbox - WP8 - BacHBerry project, co-funded by the European
% Commission in the 7th Framework Programme (Project No. FP7-613793)
% References: Koziol FEBS letters 584 (2010) 941-944; Eisinga FEBS letters 587 (2013) 677-682
%
% Susana Vinga, 2016.05.03 susanavinga@tecnico.ulisboa.pt

[n,k]=size(M);   % n -  deletions, c - columns/criteria

%sort each column descending - best (higher) on top! if some criteria is
%better when lower (for example, distance to Wild-Type, pre-process with
%symetric, e.g. M(:,4)=-M(:,4);
[Y,I]=sort(M,1,'descend');

% get indicator I - where (position) are my deletions?
[X,R]=sort(I,1,'ascend');    %gives directly the ranking R


RP = prod(R,2);

s = -log(RP)+k*log(n+1);  %test statistics

%compare with Gamma(k,1) approximation, p-value with Gamma - PG
% Koziol FEBS letters 584 (2010) 941-944

PG = 1-gamcdf(s,k,1);


% exact p-value - PE
% Eisinga FEBS letters 587 (2013) 677-682

% ---------------------------------------------------------
%--------- slower !!!-----------
% ---------------------------------------------------------

% for i=1:n %for all rows/cases
%     i
%     for j=1:RP(i) %for all values from 1 up to the RP obtained in that row
%         h(j)=piltzcount(j,k,n); %calculate probability of obtaining that RP
%     end
%     H(i)=sum(h); %sum all the probabilities
%     clear h;
%     
% end
% 
% PE=H/n^k;   %exact probabilities

% ---------------------------------------------------------
%---------other option, not to repeat values!!!-----------
% ---------------------------------------------------------

wait=1e5;  %maximum rank product allowed to be calculated (processing time restrictions)

maxR=min(max(RP),wait);
for i=1:maxR
    %i
    ha(i)=piltzcount(i,k,n);    %calculate all permutations for that rank value
end

for i=1:n
    if RP(i)<=maxR   % low rank --> calculate exact value
        H(i)=sum(ha(1:RP(i)));  % sum all probabilites up to the RP obtained
    else
        H(i)=0;    %will use gamma approximation afterwards, if necessary
    end
end

PE=H/n^k;   % exact probability

x=find(~PE);   %where are the zeros...
PE(x)=PG(x);   %...to assign gamma approximation for those


%-----TODO------bounds - to run faster
% Heskes BMC Bioinformatics 2014, 15:367
%.....................
% only available for R
%.....................
%-----TODO------bounds


% multiple testing correction - q-values - Storey (2002)

[FDR, QG] = mafdr(PG);
[FDR, QE] =mafdr(PE);

P.PG=PG;
P.PE=PE';
P.QG=QG;
P.QE=QE';