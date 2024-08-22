function [X]=tnormrnd(leftT,rightT,B,mu,tau,burnin)
% function [X]=tnormrnd(tailI,B,m,t,burnin)
%
% To sample from a truncated standard univariate normal distribution.
% Using Gibbs sampler and latent variable technique from Damien, Walker
% (2011)Sampling Truncated Normal, Beta, and Gamma Densities Journal of 
% Computational and Graphical Statistics, Vol. 10, No. 2, pp. 206-215
% 
% INPUT:   
%   leftT: indicator for left tail included in the support
%   rightT: indicator for right tail included in the support
%   B: vector of real valued truncation bounds ordered from left to right
%   (mu,tau): mean and precision parameters for the truncated Normal
%   burnin: number of iterations for the burn in period
%
% OUTPUT: 
% X: single realization of the truncated normal.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b=length(B);
B=(B-mu)*sqrt(tau); % Standardized bounds
% Calculate number of intervals in the support
p=(b+leftT+rightT)/2;
% if leftT+rightT~=1
%     p=b/2+leftT*rightT;
% else
%     p=(b-1)/2+1;
% end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ALLOCATE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P=zeros(p,1);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE PROBABILITIES FOR EACH INTERVAL %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if leftT % There is a left tail
    P(1)=normcdf(B(1)); % Left tail probability
end
if rightT % There is a right tail
    P(end)=normcdf(-B(b)); % Left tail probability
end
iB=B(1+leftT:b-rightT);
for i=1:length(iB)/2
    if iB(2*i-1)>0
        P(i+leftT)=normcdf(-iB(2*i-1))-normcdf(-iB(2*i));
    elseif iB(2*i)>0
        P(i+leftT)=1-(normcdf(-iB(2*i))+normcdf(iB(2*i-1)));
    else
        P(i+leftT)=normcdf(iB(2*i))-normcdf(iB(2*i-1));
    end
end

%% Check if the total probability of the truncation is positive%%
% Otherwise, use a suitable aproximation
if sum(P)==0
    disp('I am approximating the truncated normal')
    mb=min(abs(B));
    if leftT % There is a left tail
        P(1)=exp(-0.5*(-B(1)-mb)^2-log(-B(1)));
        %P(1)=-exp(-0.5*B(1)^2)/B(1); % Left tail probability
    end
    if rightT % There is a right tail
        P(end)=exp(-0.5*(B(b)-mb)^2-log(B(end)));
        %P(end)=exp(-0.5*B(b)^2)/B(b); % Left tail probability
    end
    for i=1:length(iB)/2 
        if iB(2*i-1)>0 % Positive interval
            P(i+leftT)=exp(-0.5*(iB(2*i-1)-mb)^2-log(iB(2*i-1)))-...
                exp(-0.5*(iB(2*i)-mb)^2-log(iB(2*i)));
        else % Negative interval. No other option since sum(P)==0
            P(i+leftT)=exp(-0.5*(-iB(2*i)-mb)^2-log(-iB(2*i)))-...
                exp(-0.5*(-iB(2*i-1)-mb)^2-log(-iB(2*i-1)));
        end
        %P(i+leftT)=exp(-0.5*iB(2*i)^2)/iB(2*i-1)-exp(-0.5*iB(2*i-1)^2)/iB(2*i-2);
    end
end
clear iB

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAMPLE FROM STANDARD TRUNCATED NORMAL%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Choose interval for sampling
i=discreternd(cumsum(P),1);

%% Initialize X an bounds for the interval
if i==1&& leftT % Sample from left tail
    type='left';
    B=B(1);
    X=B-1;
elseif i==p && rightT
    type='right';
    B=B(b);
    X=B+1;
else
    type='between';
    B=[B(2*i-(leftT+1)),B(2*i-leftT)];
    X=mean(B);
end
% Gibbs sampler
for i=1:burnin+1
   W=sqrt(-2*log(rand)+X^2);
    if strcmp(type,'right')
        l=max(B,-W);
        L=W;
    elseif strcmp(type,'left')
        l=-W;
        L=min(B,W);
    else
        l=max(B(1),-W);
        L=min(B(2),W);
    end
    X=l+(L-l)*rand;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UNSTANDARIZE%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X=X*sqrt(1/tau)+mu;