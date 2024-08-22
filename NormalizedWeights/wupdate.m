function wa=wupdate(PPType,b,Ja,n,da,ka,Da)
% function [w]=wupdate(PPType,b,Ja,n,da,ka,Da)
% Updates the mixture weights (w)
% INPUT:
% PPType: Type of Nonparametric prior may be:
%       'StickBreaking' (default) for a two parameter Stick-Breaking Prior
%       'Geometric' for Geometric Stick breaking Prior 
%       'DirichletHP' for Dirichlet process with Gamma hyperprior for the
%       mass parameter
% b: If PPType is 'StickBreaking' or 'DirichletHP', then 2 x 1 vector of 
%   hyperparameters Beta random variables associated with 
%   Stick-Breaking prior weights.
%   If PPType is 'DirichletHP' then current mass parameter of the Dirichlet process
% Ja: current number of components
% n: sample size
% da: 1 x n vector with current state of index variables
% ka: 1 x n vector with current k values
% Da: cell array of length n with current index values for latent variables
%
% OUTPUT:
% wa: vector of length Ja of updated mixture weights

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wa=zeros(1,Ja);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TWO PARAMETER STICK-BREAKING PRIOR %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(PPType,'StickBreaking')
    for j=1:Ja
        N_j=0;
        M_j=0;
        for i=1:n
            N_j=N_j+sum(Da{i}==j);
            M_j=M_j+sum(Da{i}>j);
        end
        b1j=b(1)+sum(da==j)+N_j;
        b2j=b(2)+sum(da>j)+M_j;
        vnew=betarnd(b1j,b2j);
        if j==1
            wa(1)=vnew;
        else
            wa(j)=wa(j-1)*vnew*(1-vold)/vold;
        end
        vold=vnew;
    end
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GEOMETRIC STICK-BREAKING PRIOR %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(PPType,'Geometric')
    %wa=GeomPupdate_weights(b(1),b(2),Ja,n,da,ka,Da);
    b1hat=b(1)+n;
    b2hat=b(2)-n;
    for i=1:n
        b1hat=b1hat+ka(i);
        b2hat=b2hat+da(i)-ka(i);
        if ka(i)>0
            for j=1:ka(i)
                b2hat=b2hat+Da{i}(j);
            end
        end
    end
    lambda=betarnd(b1hat,b2hat);
    wa=lambda*(1-lambda).^((0:Ja-1)');
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DIRICHLET PRIOR WITH GAMMA HYPERPRIOR FOR THE MASS PARAMETER %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(PPType,'DirichletHP')
    for j=1:Ja
        N_j=0;
        M_j=0;
        for i=1:n
            N_j=N_j+sum(Da{i}==j);
            M_j=M_j+sum(Da{i}>j);
        end
        b1j=1+sum(da==j)+N_j;
        b2j=b+sum(da>j)+M_j;
        vnew=betarnd(b1j,b2j);
        if j==1
            wa(1)=vnew;
        else
            wa(j)=wa(j-1)*vnew*(1-vold)/vold;
        end
        vold=vnew;
    end
end
