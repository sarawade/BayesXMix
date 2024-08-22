function [mua]=muupdate(x,n,q,p,da,Ja,taua,ka,Da,U,u,mu0,c,nitn)
% function [mua]=muupdate(x,n,q,p,da,Ja,taua,ka,Da,U,u,mu0,c,nitn)
% Updates the current state of the mean parameters for continuous 
% covariates in NPRegNW (normal Kernel and normal base measure)
% INPUT:
% x: n x p matrix of continuous covariates 
% n: sample size
% q: number of discrete covariates 
% p: number of continuous covariates 
% da: 1 x n vector with current state of index variables
% Ja: current number of components
% taua: p x 1 vector with the current state of the common component
%   precision
% ka: 1 x n vector with current k values
% Da: cell array of length n with current index values for latent variables
% (U,u): cell arrays of length n with current uniform latent variables and
%   associated indicators
% [mu0,c]: hyperparameters for the Normal prior for the (mu_j|tau)
% nitn: number of iterations for truncated normal sampling
%
% OUTPUT:
% mua: p x Ja matrix of updated mean parameters

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ALLOCATION %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mua=zeros(p,Ja);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAMPLING FOR EACH COMPONENT %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:Ja
    %% Calculate posterior parameters
	n_j=sum(da==j); % number of observations associated to component j
    if n_j>0
        c_hat=c+n_j; % updated c
        mu0_hat=1./c_hat.*(c.*mu0+ sum(x(da==j,:),1)'); % updated mu0
        t_hat=c_hat.*taua; % posterior precision
    else % No observations are associated to component j, so use prior
        mu0_hat=mu0;
        t_hat=c.*taua;
    end
    %% Sample each component of truncated multivariate normal
    for h=1:p
        %% Find the truncation for mu(h,j)
        [T,B]=mutruncation(j,h,x(:,h),n,q,taua(h),ka,Da,U,u);
        %% Sample from posterior
        if isempty(B) % No truncation so sample from a normal
            mua(h,j)=mu0_hat(h)+randn*sqrt(1/t_hat(h));
        else % Sample from a truncated normal
            mua(h,j)=tnormrnd(T,T,B,mu0_hat(h),t_hat(h),nitn-1);
        end
    end
end