
%   The exact probability distribution of the rank product statistics for 
%   replicated experiments, by Eisinga, Breitling & Heskes, FEBS Letters, 
%   January 2013


function h = piltzcount(r,k,n)

%PILTZCOUNT Number of ways to construct a rank product
%   H = piltzcount(R,K,N) returns the number of ways rank product R can be
%   constructed from K experiments if N is the number of genes (and thus
%   the maximum rank in a single experiment)
%
%   Implements MYFACTOR, MYDIVISOR, MYNCHOOSEK, and MYMULTCOEF


if nargin < 3,
    n = r;
end

% Catch all trivial cases

if k < 1,
    h = 0;
    return;
elseif r == 1,
    h = 1;
    return;
elseif r > n^k,
    h = 0;
    return;
elseif k == 1,
    h = 1;
    return;
end

% Compute prime factorization

[f,m] = myfactor(r);
nprimes = length(f);

% First term: number of ordered k-tuples such that their product equals r

h = 1;
for t=1:nprimes,
    h = h*mynchoosek(m(t)+k-1,k-1);
end

% Now consider additional terms

if r > n,
    % Get all divisors of r larger than n
    
    d = mydivisor(f,m);
    bigones = d(d > n);
    nbig = length(bigones);
    
    % Start constructing possible combinations of divisors
    
    divset = cell(1,k);
    bbb = cell(1,k);
    divset{1} = 1;
    bbb{1} = zeros(1,nbig);
    smax = ceil(log(r)/log(n)) - 1;
                     % maximum number of divisors that can be divided out
    for s=1:smax,
        newdivset = divset{s}(:)*bigones(:)';
                     % any combination of s+1 divisors consists of an
                     % allowed combination of s divisors plus another big
                     % divisor
        [i,j] = find(r - newdivset.*floor(r./newdivset) < 0.5);
                     % checks whether the remainder is an integer
        nnew = length(i);
        if nnew,
            divset{s+1} = zeros(1,nnew);
            bbb{s+1} = zeros(nnew,nbig);
            for t=1:nnew,
                divset{s+1}(t) = divset{s}(i(t))*bigones(j(t));
                bbb{s+1}(t,:) = bbb{s}(i(t),:);
                bbb{s+1}(t,j(t)) = bbb{s+1}(t,j(t))+1;
                     % new combination t at level s+1 is combination i(t)
                     % at level s plus divisor j(t)
            end
            [bbb{s+1},iii] = unique(bbb{s+1},'rows');
                     % some combinations can be constructed through
                     % different routes: keep them only once
            divset{s+1} = divset{s+1}(iii);
            for t=1:length(iii),
                hextra = piltzcount(r/divset{s+1}(t),k-s,r);
                     % number of ordered k-s-tuples such that their product
                     % equals r/divset
                h = h + (-1)^s * mynchoosek(k,s) * mymultcoef(s,bbb{s+1}(t,:)) * ...
                    hextra;
            end
        else
            s = smax;
                     % if there are no combinations left at level s, there
                     % can be none at higher levels
        end
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f,m] = myfactor(n)

%MYFACTOR Prime factors and multiplicities
%    F = MYFACTOR(N) returns the primes in the prime factorization of N
%    [F,M] = MYFACTOR(N) also returns the corresponding multiplicities
%    Example: [f,m] = myfactor(12); yields f = [2,3] and m = [2,1]

% Catch trivial case

if n == 1,
    f = 1;
    m = 0;
    return;
end

% Use Matlab's FACTOR to get the prime factorization

pf = factor(n);
f = unique(pf);

% Compute multiplicities

if nargout > 1,
    nprimes = length(f);
    m = zeros(1,nprimes);
    for i=1:nprimes,
        m(i) = sum(pf == f(i));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function d = mydivisor(f,m)

%MYDIVISOR List of integer divisors of a number
%    D = MYDIVISOR(F,M) returns a row vector D of all distinct divisors
%    of a positive integer making use of the primes F and corresponding
%    multiplicities M in its prime factorization

% Each divisor can be written as a product of prime factors

nprimes = length(f);
d = f(1).^(0:1:m(1))';
for i = 2:nprimes,
    d = d*(f(i).^(0:1:m(i)));
    d = d(:);
end
d = sort(d)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = mynchoosek(n,k)

%MYNCHOOSEK Efficient calculation of binomial coefficient

y = exp(gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1));
y = round(y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = mymultcoef(n,x)

%MYMULTCOEF Efficient calculation of multinomial coefficient
%    X should be a vector of integers adding up to N

y = exp(gammaln(n+1)-sum(gammaln(x+1)));
y = round(y);


