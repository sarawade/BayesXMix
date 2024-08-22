function [da,Da,sigmaa,betaa,rhoa,mua,wa]=move2(n,q,p,sigmaa,betaa,rhoa,mua,wa,da,Da)
% [da,Da,sigmaa,betaa,rhoa,mua,wa]=move2(n,q,p,sigmaa,betaa,rhoa,mua,wa,da,Da);
% Proposes first label-switching move
% INPUT:
% n: sample size
% q: number of discrete covariates 
% p: number of continuous covariates 
% betaa: q+p+1 x Ja matrix of current regression coefficients
% sigmaa: 1 x Ja vector of current regression variance parameters
% rhoa: cell array of length q, where cell h is a matrix of size Gh x Ja of
%   current discrete covariate parameters
% mua: p x Ja matrix of current mean parameters
% wa: vector of length Ja of current mixture weights
% da: 1 x n vector with current state of index variables
% Da: cell array of length n with current index values for latent variables

% OUTPUT:
% betaa: q+p+1 x Ja matrix of updated regression coefficients
% sigmaa: 1 x Ja vector of updated regression variance parameters
% rhoa: cell array of length q, where cell h is a matrix of size Gh x Ja of
%   updated discrete covariate parameters
% mua: p x Ja matrix of updated mean parameters
% wa: vector of length Ja of updated mixture weights
% da: 1 x n vector with updated state of index variables
% Da: cell array of length n with updated index values for latent variables

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAMPLE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Choose one component among 1 to J_occup-1 with uniform probability
%% where J_occup is the index of the largest occupied component

% Compute the largest occupied component
J_occup=max(da);
for i=1:n
    J_occup=max([J_occup, Da{i}]);
end

% Uniformly choose one
j=discreternd(1:(J_occup-1),1);

%% Swap labels, parameters, and weights associated to components j and j+1
% Labels
dp=da;
Dp=Da;
dp(da==j)=j+1;
dp(da==(j+1))=j;
for i=1:n
    Dp{i}(Dp{i}==j)=j+1;
    Dp{i}(Dp{i}==(j+1))=j;
end
% Parameters
betap=betaa;
sigmap=sigmaa;
betap(:,j)=betaa(:,j+1);
betap(:,j+1)=betaa(:,j);
sigmap(j)=sigmaa(j+1);
sigmap(j+1)=sigmaa(j);
if q>0&&p>0
    rhop=rhoa;
    mup=mua;
    mup(:,j)=mua(:,j+1);
    mup(:,j+1)=mua(:,j);
    for hq=1:q
        rhop{q}(:,j)=rhoa{q}(:,j+1);   
        rhop{q}(:,j+1)=rhoa{q}(:,j);
    end
elseif q>0
    rhop=rhoa;
    for hq=1:q
        rhop{q}(:,j)=rhoa{q}(:,j+1);   
        rhop{q}(:,j+1)=rhoa{q}(:,j);
    end
else
    mup=mua;
    mup(:,j)=mua(:,j+1);
    mup(:,j+1)=mua(:,j);
end
% Weights
wp=wa;
Vjp1=wa(j+1)/(1-sum(wa(1:j)));
Vj=wa(j);
if j>1
    Vj=Vj/(1-sum(wa(1:(j-1))));      
end    
wp(j)=Vjp1;
if j>1
    Vjm1=wa(j-1);
    if j>2
        Vjm1=Vjm1/(1-sum(wa(1:(j-2))));
    end
    wp(j)=wp(j)*wa(j-1)*(1-Vjm1)/Vjm1;
end
wp(j+1)=wp(j)*Vj*(1-Vjp1)/Vjp1;

%% Compute acceptance probability
nj=sum(da==j);
njp1=sum(da==(j+1));

Nj=0;
Njp1=0;
for i=1:n
    Nj=Nj+sum(Da{i}==j);
    Njp1=Njp1+sum(Da{i}==(j+1));
end

accep=min(1,(1-Vjp1)^(nj+Nj)/((1-Vj)^(njp1+Njp1)));
if rand<accep
    da=dp;
    Da=Dp;
    wa=wp;
    betaa=betap;
    sigmaa=sigmap;
    if q>0&&p>0
    rhoa=rhop;
    mua=mup;
    elseif q>0
        rhoa=rhop;
    else
        mua=mup;
    end
end
