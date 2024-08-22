function p=tdirichletrnd(a,m,M,p,burnin)
% function [p]=tdirichletrnd(a,m,M,p,burnin)
% Generate a realization from a G-dimensional truncated Dirichlet distribution with parameter a
%
% INPUT
%   a: G x 1 vector of parameters for the Dirichlet distribution. 
%   (m,M): G x 1 vectors of lower and upper bounds respectively
%   p: initial value for Gibbs sampling. MUST satisfy truncation (and add up
%       to 1)
%   burnin: number of iterations for the burn in period
% 
% OUTPUT
%   p: G x 1 realization from the truncated Dirichlet distribution.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G=size(a,1);
Z=p*gamrnd(sum(a),1); % Latent gamma variables
B=zeros(1,2); % Vector of bounds for the full conditional distributions

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VALIDATE STARTING POINT %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for g=1:G
    if p(g)>M(g) || p(g)<m(g)
        disp('Warning: initial point is outside truncation. Error may occur')
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAMPLE%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Latent variables
for i=1:burnin+1
    for g=1:G
        % Calculate full conditional lower bound
        B(1)=max((1-M(1:G~=g))./M(1:G~=g).*Z(1:G~=g)-(ones(G-1)-eye(G-1))*Z(1:G~=g));
        B(1)=max(B(1),m(g)/(1-m(g))*sum(Z(1:G~=g)));
        if max(m(1:G~=g))==0 && M(g)==1 % No upper bound
            if B(1)==0 % No lower bound
                Z(g)=gamrnd(a(g),1);
            else % Lower bound only so sample from right tail
                Z(g)=tgamrnd('right',B(1),a(g),1,Z(g),0);
            end
        else
            % Calculate full conditional upper bound
            B(2)=min((1-m(1:G~=g))./m(1:G~=g).*Z(1:G~=g)-(ones(G-1)-eye(G-1))*Z(1:G~=g));
            B(2)=min(B(2),M(g)/(1-M(g))*sum(Z(1:G~=g)));
            if B(1)==0 % Upper bound only so sample from left tail
                Z(g)=tgamrnd('left',B(2),a(g),1,Z(g),0);
            else % Both bounds so sample from bounded interval
                Z(g)=tgamrnd('between',B,a(g),1,Z(g),0);
            end
        end
    end
end

%% calculate p by standardizing latent variables
p=Z/sum(Z);