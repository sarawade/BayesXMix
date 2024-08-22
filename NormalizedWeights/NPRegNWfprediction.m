function [Y_fpred,lpred,upred]=NPRegNWfprediction(Y_grid,x_new,q,p,w,theta,J,varargin)
%function [Y_pred,lpred,upred]=NPRegNWfprediction(Y_grid,x_new,q,p,w,theta,J,credprob)
%
% Compute estimated predictive density
% INPUT:
%   Y_grid: column vector of Y values at which predictive density will be
%       estimated
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
%   Y_fpred: n_new x length(Y_grid) matrix of predictive density estimates.
%   lpred, upred: n_new x length(Y_grid) matrices of lower and upper bounds 
%                 for the pointwise credible intervals. Only returned if 
%                 cred_prob argument is specified

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE AND ALLOCATE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ny=length(Y_grid); 
n_new=size(x_new,1); % Number of new points for prediction
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
Y_fpred=zeros(n_new,ny);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE PREDICTIVE DENSITY ESTIMATES %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if ~isempty(varargin)
    for s=1:S
        % Compute covariate dependent weight matrix
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
        for i=1:n_new
            Y_fpred(i,:)=Y_fpred(i,:)+sum(repmat(wx(i,:)',1,ny).*normpdf(repmat(Y_grid',J(s),1),repmat(betax(i,:)',1,ny),repmat(sqrt(sigma{s})',1,ny)),1);
        end
    end
    % MCMC predictive estimate
    Y_fpred=Y_fpred/S; 
%end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE CREDIBLE INTERVALS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(varargin) % Optional argument is used, so intervals are calculated
    lpred=zeros(n_new,ny);
    upred=zeros(n_new,ny);
    for i=1:n_new
        fpredS=zeros(S,ny);
        for s=1:S
            % Compute covariate dependent weight matrix
            lg1=zeros(1,J(s));
            lg2=zeros(1,J(s));
            if q>0
                for h=1:q
                    T=size(rho{s}{h},1);
                    lg1=lg1+((ones(1,T)*x_new(i,h))==(0:T-1))*log(rho{s}{h});
                end
            end
            if p>0
                for h=1:p
                    lg2=lg2-0.5*tau(s,h)*((ones(1,J(s))*x_new(i,q+h))-mu{s}(h,:)).^2;
                end
            end
            wx=w{s}.*exp(lg2+lg1); % Weights *normal kernel for x_new in component j
            wx=wx/(sum(wx)); % Normalize
            % Compute prediction matrix
            betax=[1,x_new(i,:)]*beta{s}; % Predictive mean for x_new in component j
            fpredS(s,:)=sum(repmat(wx',1,ny).*normpdf(repmat(Y_grid',J(s),1),repmat(betax',1,ny),repmat(sqrt(sigma{s})',1,ny)),1);
        end       
        % MCMC predictive estimate
        Y_fpred(i,:)=sum(fpredS)/S; 
        cred_prob=varargin{1};
        fpredS_sort=sort(fpredS); 
        lpred(i,:)=fpredS_sort(uint16(S*(1-cred_prob)/2),:); % Lower bound
        upred(i,:)=fpredS_sort(uint16(S*(cred_prob+(1-cred_prob)/2)),:); % Upper bound
    end
end
