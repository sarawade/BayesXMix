function [x]=discreternd(p,n)
% function [x]=discreternd(p,n)
% To sample n iid variables from a discrete distribution with P(X<=i) proportional to p(i)
%
% INPUT:
%   p: m x 1 vector of CUMULATIVE probabilities (not necessarily  standardized)
%   n: sample size
% 
% OUTPUT:
%   x: n x 1 vector of realizations

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m=length(p);
x=ones(n,1);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAMPLE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:n
    u=rand;
    while u>p(x(i))/p(m)
        x(i)=x(i)+1;
    end
end