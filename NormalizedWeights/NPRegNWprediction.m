function [Y_pred,predS,Y_pred_CI]=NPRegNWprediction(x_new,q,p,w,theta,J,varargin)
%function [Y_pred,Y_pred_CI]=NPRegNWprediction(x_new,q,p,w,theta,J,cred_prob)
%
% Compute prediction with pointwise credible intervals
% INPUT:
%   x_new: n_new x q+p matrix of new covariates
%   q: number of discrete covariates
%   p: number of continuous covariates
%   w: cell array of length S of posterior samples of mixture weights
%   theta: cell array of length S of posterior samples of mixture
%       parameters
%   J: S x 1 vector of posterior samples of number of mixture components
%       with positive weights
%
% OPTIONAL ARGUMENTS:
%   cred_prob: Probability for pointwise credible intervals. 
%
% OUTPUT:
%   Y_pred: n_new x 1 vector of predictive mean of Y|X_new,data
%   Y_pred_CI: n_new x 2 matrix of lower and upper bounds for the pointwise
%       credible intervals. Only returned if cred_prob argument is specified

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE AND ALLOCATE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_new=size(x_new,1); % Number of new points for prediction
S=length(J); % MCMC posterior sample size
% Model parameters posterior sample
beta=theta{1};
if q==0 % Only continuous variables
    mu=theta{3};
    tau=theta{4};
else % Discrete variables
    rho=theta{3};
    if p>0 % Continuous variables
        mu=theta{4};
        tau=theta{5};
    end
end
predS=zeros(S,n_new);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE PREDICTIVE ESTIMATES %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for s=1:S
    % Compute covariate dependent weights
    lg1=zeros(n_new,J(s));
    lg2=zeros(n_new,J(s));
    if q>0
        for h=1:q
            T=size(rho{s}{h},1);
            lg1=lg1+(kron(ones(1,T),x_new(:,h))==kron(ones(n_new,1),(0:T-1)))*log(rho{s}{h});
        end
    end
    if p>0
        for h=1:p
            lg2=lg2-0.5*tau(s,h)*(kron(ones(1,J(s)),x_new(:,q+h))-kron(ones(n_new,1),mu{s}(h,:))).^2;
        end
    end
    wx=kron(ones(n_new,1),w{s}).*exp(lg2+lg1); % Weights *normal kernel for x_new in component j
    wx=wx./kron(ones(1,J(s)),sum(wx,2)); % Normalize
    % Compute prediction matrix
    betax=[ones(n_new,1),x_new]*beta{s}; % Predictive mean for x_new in component j
    predS(s,:)=sum(wx.*betax,2)'; % Average the prediction over the mixture components
end
% MCMC predictive estimate
Y_pred=(sum(predS,1)/S)'; 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE CREDIBLE INTERVALS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(varargin) % Optional argument is used, so intervals are calculated
    cred_prob=varargin{1};
    predS_sort=sort(predS);
    m_pred=predS_sort(uint16(S/2),:); % Median
    l_pred=predS_sort(uint16(S*(1-cred_prob)/2),:); % Lower bound
    u_pred=predS_sort(uint16(S*(cred_prob+(1-cred_prob)/2)),:); % Upper bound
    Y_pred_CI=[l_pred',u_pred',m_pred'];
end


