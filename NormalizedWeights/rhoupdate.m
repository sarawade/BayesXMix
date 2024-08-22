function [rhoa]=rhoupdate(x,n,q,da,ka,Da,U,u,Ja,rhoa,gamma,nitg)
% function [rhoa]=rhoupdate(x,n,q,da,ka,Da,U,u,Ja,rhoa,gamma,nitg)
% Updates the current state of the parameters for discrete 
% covariates in NPRegNW (Multinomial kernel and Dirichlet base measure)
% INPUT:
% x: n x q matrix of discrete covariates 
% n: sample size
% q: number of discrete covariates
% da: 1 x n vector with current state of index variables
% ka: 1 x n vector with current k values
% Da: cell array of length n with current index values for latent variables
% (U,u): cell arrays of length n with current uniform latent variables and
%   associated indicators
% Ja: current number of components
% rhoa: cell array of length q, where cell h is a matrix of size Gh x Ja of
%   current discrete covariate parameters
% gamma: hyperparameters for the Dirichlet prior for the (rho_j)
% nitg: number of iterations for truncated gamma sampling
%
% OUTPUT:
% rhoa: cell array of length q, where cell h is a matrix of size Gh x Ja of
%   updated discrete covariate parameters

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAMPLING FOR EACH COVARIATE AND EACH COMPONENT  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for h=1:q
    Gh=size(gamma{h},1);
    for j=1:Ja    
        %% Calculate posterior parameters 
        gammah_hat=gamma{h}+sum(kron(ones(1,Gh),x(da==j,h))==kron(0*x(da==j,h)+1,0:Gh-1),1)';
        %% Find the truncation for rho(h,j) (Gh dimensional columns)
        [lbound,ubound]=rhotruncation(j,h,Gh,x(:,h),n,ka,Da,U,u);
        %% Sample from the posterior
        if sum(ubound)==Gh && sum(lbound)==0 % No truncation so sample from Dirichlet
            rhoa{h}(:,j)=dirichletrnd(gammah_hat,1);
        else % Sample from truncated Dirichlet
            rhoa{h}(:,j)=tdirichletrnd(gammah_hat,lbound,ubound,rhoa{h}(:,j),nitg-1);
        end
    end
end
