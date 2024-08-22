function [lb,ub,ubI]=tautruncation(x,n,q,p,mua,ka,Da,U,u)
% function [lb,ub,ubI]=tautruncation(x,n,q,p,mua,ka,Da,U,u)
% Calculates truncation for posterior sampling of the precision parameter 
% for continuous covariates in NPRegNW (normal Kernel and Gamma prior)
% INPUT:
% x: n x p matrix of continuous covariates 
% n: sample size
% q: number of discrete covariates 
% p: number of continuous covariates 
% mua: p x Ja matrix of current mean parameters
% ka: 1 x n vector with current k values
% Da: cell array of length n with current index values for latent variables
% (U,u): cell arrays of length n with current uniform latent variables and
%   associated indicators
%
% OUTPUT:
% lb: p x 1 vector of lower bounds
% ub: p x 1 vector of upper bounds (=0 if no upper bound)
% ubI: p x 1 vector of indicators of upper bound presence

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lb=zeros(p,1);
ub=zeros(p,1);
ubI=zeros(p,1);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIND TRUNCATION INTERVALS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:n
    if ka(i)>0
        for l=1:ka(i)
            vil=((-2*log(U{i}(l,q+1:end)))./((x(i,:)-mua(:,Da{i}(l))').^2))';
            aux=u{i}(l,q+1:end);
            lb(aux==0)=max(lb(aux==0), vil(aux==0));
            ub(aux==1)=ubI(aux==1).*min(ub(aux==1),vil(aux==1))+(ubI(aux==1)==0).*vil(aux==1);
            ubI(aux==1)=1;
        end
    end
end
