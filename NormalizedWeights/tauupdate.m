function [taua]=tauupdate(x,n,q,p,da,Ja,mua,taua,ka,Da,U,u,mu0,c,a1,a2,nitg)
% function [taua]=tauupdate(x,n,q,p,da,Ja,mua,taua,ka,Da,U,u,mu0,c,a1,a2,nitg)
% Updates the current state of the precision parameter for continuous 
% covariates in NPRegNW (normal Kernel and Gamma prior)
% INPUT:
% x: n x p matrix of continuous covariates 
% n: sample size
% q: number of discrete covariates 
% p: number of continuous covariates 
% da: 1 x n vector with current state of index variables
% Ja: current number of components
% mua: p x Ja matrix of current mean parameters
% taua: p x 1 vector with the current state of the common component
%   precision
% ka: 1 x n vector with current k values
% Da: cell array of length n with current index values for latent variables
% (U,u): cell arrays of length n with current uniform latent variables and
%   associated indicators
% [mu0,c,a1,a2](if p>0): hyperparameters for the Normal-Gamma prior for
%   the (mu_j,tau). 
% nitg: number of iterations for truncated Gamma sampling
%
% OUTPUT:
% taua: p x 1 vector with the updated state of the common component
%   precision

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE POSTERIOR PARAMETERS AND TRUNCATION %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Posterior parameters
a2_hat=a2+.5*(sum(x.^2,1)');
a1_hat=a1+Ja/2;
for j=1:Ja
	n_j=sum(da==j); % number of observations associated to component j
    if n_j>0
        c_hat=c+n_j; % updated c
        mu0_hat=1./c_hat.*(c.*mu0+ sum(x(da==j,:),1)'); % updated mu0
        a2_hat=a2_hat+.5*(c.*mu0.^2-c_hat.*mu0_hat.^2 )+.5*(c_hat.*(mua(:,j)-mu0_hat).^2);
    end
end

%% Find the truncation %%
[lbound,ubound,uboundIndicator]=tautruncation(x,n,q,p,mua,ka,Da,U,u);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAMPLING FOR EACH COVARIATE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for h=1:p
    if lbound(h)==0 && uboundIndicator(h)==0
        %% No truncation, sample from a gamma
        taua(h)=gamrnd(a1_hat(h),1/a2_hat(h));
    elseif lbound(h)==0
        %% Right Truncation, sample from left tail of a gamma
        taua(h)=tgamrnd('left',ubound(h),a1_hat(h),a2_hat(h),taua(h),nitg-1);
    elseif uboundIndicator(h)==0
        %% Left Truncation, sample from right tail of a gamma
        taua(h)=tgamrnd('right',lbound(h),a1_hat(h),a2_hat(h),taua(h),nitg-1);
    else
        %% Both Truncations, sample in an interval
        taua(h)=tgamrnd('between',[lbound(h),ubound(h)],a1_hat(h),a2_hat(h),taua(h),nitg-1);
    end
end
