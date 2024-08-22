function [Jd,JD,Ja,rhoa,mua]=Jupdate(phi,n,q,p,da,ka,Da,Ja,rhoa,mua,taua,gamma,mu0,c)
% function [Jd,JD,Ja,rhoa,mua]=Jupdate(phi,n,q,p,da,ka,Da,Ja,rhoa,mua,taua,gamma,mu0,c)
% Updates the number of components that are needed for MCMC update
% INPUT:
% phi: Slice sampling constant
% n: sample size
% q: number of discrete covariates 
% p: number of continuous covariates 
% da: 1 x n vector with current state of index variables
% ka: 1 x n vector with current k values
% Da: cell array of length n with current index values for latent variables
% Ja: current number of components
% rhoa: cell array of length q, where cell h is a matrix of size Gh x Ja of
%   current discrete covariate parameters
% mua: p x Ja matrix of current mean parameters
% taua: p x 1 vector with the current state of the common component
%   precision
% gamma (if q>0): hyperparameters discrete covariate parameters
% mu0,c: hyperparameters for the continuous covariate parameters
%
% OUTPUT:
% Jd: updated vector of number of components needed to sample the d
% JD: updated cell array of number of components needed to sample the D
% Ja: updated total number of components
% (rhoa,mua): with updated size corresponding to the new Ja

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jd=zeros(n,1);
JD=cell(n,1);
J=0;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAMPLE J'S %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:n
    Jd(i)=da(i)+fix(-log(rand)/phi); 
    J=max(J,Jd(i)); 
    JD{i}=zeros(ka(i),1);
    if ka(i)>0
        for j=1:ka(i)
            JD{i}(j)=Da{i}(j)+fix(-log(rand)/phi); 
            J=max(J,JD{i}(j));
        end
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAMPLE NEW COMPONENT PARAMETERS IF REQUIRED %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if J>Ja
    if q>0
        for h=1:q
            rhoa{h}(:,Ja+1:J)=dirichletrnd(gamma{h},J-Ja);
        end
    end
    if p>0
        mua(:,Ja+1:J)=randn(p,J-Ja)./kron(ones(1,J-Ja),(taua.*c).^.5)+kron(ones(1,J-Ja),mu0);
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REMOVE COMPONENTS IF REQUIRED %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if J<Ja
    if q>0
        for h=1:q
            rhoa{h}=rhoa{h}(:,1:J);
        end
    end
    if p>0
        mua=mua(:,1:J);
    end
end

Ja=J;