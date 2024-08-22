function [p]=dirichletrnd(a,n)
% function [p]=dirichletrnd(a,n)
% Generate n realizations from a G-dimensional Dirichlet distribution with parameter a
%
% INPUT
%   a: G x 1 vector of parameters for the Dirichlet distribution. 
%   n: sample size
% 
% OUTPUT
%   p: G x n matrix of realizations from the Dirichlet distribution.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G=length(a);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAMPLE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z=gamrnd(kron(ones(1,n),a),1);
p=Z./kron(ones(G,1),sum(Z));