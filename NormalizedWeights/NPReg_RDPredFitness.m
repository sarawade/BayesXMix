function [Discrepancy_Y,Y_rep,Discrepancy_Y_rep,estimated_p]=NPReg_RDPredFitness(Y,x,q,p,w,theta,J)
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
n=size(Y,1); % Sample size
S=length(J); % MCMC posterior sample size
% Model parameters posterior sample
beta=theta{1};
sigma=theta{2};
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
Y_rep=zeros(n,S); % Replicate samples
Discrepancy_Y=zeros(S,1);
Discrepancy_Y_rep=zeros(S,1);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE PREDICTIVE ESTIMATES %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for s=1:S
    % Compute covariate dependent weights
    lg1=zeros(n,J(s));
    lg2=zeros(n,J(s));
    if q>0
        for h=1:q
            T=size(rho{s}{h},1);
            lg1=lg1+(kron(ones(1,T),x(:,h))==kron(ones(n,1),(0:T-1)))*log(rho{s}{h});
        end
    end
    if p>0
        for h=1:p
            lg2=lg2-0.5*tau(s,h)*(kron(ones(1,J(s)),x(:,q+h))-kron(ones(n,1),mu{s}(h,:))).^2;
        end
    end
    wx=kron(ones(n,1),w{s}).*exp(lg2+lg1); % Weights *normal kernel for x in component j
    wx=wx./kron(ones(1,J(s)),sum(wx,2)); % Normalize
    % Compute estimated means and variances
    betax=[ones(n,1),x]*beta{s}; % Predictive mean for x in component j
    expectation=sum(wx.*betax,2); % Average the prediction over the mixture components
    variance=wx*sigma{s}'+sum(wx.*betax.^2,2)-expectation.^2;
    % Compute observed discrepancy
    Discrepancy_Y(s)=sum((Y-expectation).^2./variance);
    % Sample replicate from the model with estimated parameters
    for i=1:n
        j=discreternd(cumsum(wx(i,:)),1);
        Y_rep(i,s)=normrnd(betax(i,j),sqrt(sigma{s}(j)));
    end
    % Computed realized discrepancy
    Discrepancy_Y_rep(s)=sum((Y_rep(:,s)-expectation).^2./variance);
end
% Calculate the posterior predictive p value
estimated_p=sum(Discrepancy_Y_rep>Discrepancy_Y)/S;



