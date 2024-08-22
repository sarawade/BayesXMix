function [lbound,ubound]=rhotruncation(j,h,Gh,x,n,ka,Da,U,u)
% function [m,M]=rhotruncation(j,h,Gh,x,n,ka,Da,U,u)
% Finds the truncation for rho{h}(:,j)
% INPUT:
% INPUT:
% j: current component
% h: current continuous covariate coordinate
% Gh: number of categories for discrete covariate h
% x: n x 1 vector of discrete covariate h  
% n: sample size
% ka: 1 x n vector with current k values
% Da: cell array of length n with current index values for latent variables
% (U,u): cell arrays of length n with current uniform latent variables and
%   associated indicators
%
% OUTPUT:
% (lbound, ubound): Gh dimensional vectors of lower and upper bounds

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ubound=ones(Gh,1);
lbound=zeros(Gh,1);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIND TRUNCATION %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:n
    if ka(i)>0
        for l=1:ka(i)
            if Da{i}(l)==j
                bilh=U{i}(l,h)*(u{i}(l,h)==1)*(x(i)==(0:Gh-1)');
                Bilh=U{i}(l,h).^((x(i)==(0:Gh-1)').*(u{i}(l,h)==0));
                lbound=max(lbound,bilh);
                ubound=min(ubound,Bilh);
            end
        end
    end
end