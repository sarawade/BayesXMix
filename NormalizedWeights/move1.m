function [da,Da,sigmaa,betaa,rhoa,mua]=move1(n,q,p,sigmaa,betaa,rhoa,mua,wa,da,Da)
% function [da,Da,sigmaa,betaa,rhoa,mua]=move1(n,q,p,sigmaa,betaa,rhoa,mua,wa,da,Da)
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
% da: 1 x n vector with updated state of index variables
% Da: cell array of length n with updated index values for latent variables

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAMPLE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Choose two non-empty components to swap with uniform probability
% Find non-empty components
uniqued=unique(da);
for i=1:n
    uniqued=unique([uniqued, Da{i}]);
end

% Uniformly choose two
ij=discreternd(1:length(uniqued),1);
j=uniqued(ij);
ih=discreternd(1:(length(uniqued)-1),1);
ih=ih+(ih>=ij);
h=uniqued(ih);

%% Swap labels and parameters associated to componenets j and h
% Labels
dp=da;
Dp=Da;
dp(da==j)=h;
dp(da==h)=j;
for i=1:n
    Dp{i}(Dp{i}==j)=h;
    Dp{i}(Dp{i}==h)=j;
end
% Parameters
betap=betaa;
sigmap=sigmaa;
betap(:,j)=betaa(:,h);
betap(:,h)=betaa(:,j);
sigmap(j)=sigmaa(h);
sigmap(h)=sigmaa(j);
if q>0&&p>0
    rhop=rhoa;
    mup=mua;
    mup(:,j)=mua(:,h);
    mup(:,h)=mua(:,j);
    for hq=1:q
        rhop{q}(:,j)=rhoa{q}(:,h);   
        rhop{q}(:,h)=rhoa{q}(:,j);
    end
elseif q>0
    rhop=rhoa;
    for hq=1:q
        rhop{q}(:,j)=rhoa{q}(:,h);   
        rhop{q}(:,h)=rhoa{q}(:,j);
    end
else
    mup=mua;
    mup(:,j)=mua(:,h);
    mup(:,h)=mua(:,j);
end

%% Compute acceptance probability
nj=sum(da==j);
nh=sum(da==h);

Nj=0;
Nh=0;
for i=1:n
    Nj=Nj+sum(Da{i}==j);
    Nh=Nh+sum(Da{i}==h);
end

accep=min(1,(wa(j)/wa(h))^(nh+Nh-nj-Nj));
if rand<accep
    da=dp;
    Da=Dp;
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
