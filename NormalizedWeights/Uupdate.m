function [U,u]=Uupdate(x,n,q,p,rhoa,taua,mua,ka,Da)
% function [U,u]=Uupdate(x,n,q,p,rhoa,taua,mua,ka,Da)
% Updates the uniform latent variables and associated indicators
% INPUT:
% x: n x q+p matrix of covariates 
% n: sample size
% q: number of discrete covariates 
% p: number of continuous covariates 
% rhoa: cell array of length q, where cell h is a matrix of size Gh x Ja of
%   current discrete covariate parameters
% taua: p x 1 vector with the current state of the common component
%   precision
% mua: p x Ja matrix of current mean parameters
% ka: 1 x n vector with current k values
% Da: cell array of length n with current index values for latent variables
%
% OUTPUT:
% (U,u): cell arrays of length n with updated uniform latent variables and
%   associated indicators

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ALLOCATION %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U=cell(n,1);
u=cell(n,1);
for i=1:n 
    % U{i}, u{i} are k(i) x q+p matrices
    if ka(i)>0 % U and u are not empty
        for l=1:ka(i)
            % Sample each p+q dimensional vector of uniform variables %%
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % CALCULATE SAMPLING SUPPORT %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            lg=zeros(q+p,1); 
            if q>0 % Discrete variable contribution
                for h=1:q
                    Gh=size(rhoa{h},1);
                    lg(h)=sum(log(rhoa{h}(:,Da{i}(l))).*((x(i,h)*ones(Gh,1))==((0:Gh-1)')));
                end
            end
            if p>0 % Continuous variable contribution
                lg(q+1:end)=-0.5*taua.*(x(i,q+1:end)'-mua(:,Da{i}(l))).^2;
            end
            
            % Create matrix of all possible combinations of lg(h) and
            % log(1-exp(lg(h)) to calculate probabilities for each of the
            % 2^(q+p)-1 sampling regions
            lP_mat=[log(1-exp(lg(1))); lg(1)]; % Matrix of log lengths 
            if (q+p)>1 % More than one covariate
                for h=2:(q+p)
                    lP_mat=[kron(ones(2,1),lP_mat), kron([log(1-exp(lg(h))); lg(h)],ones(2^(h-1),1))];
                end
            end
            P=exp(sum(lP_mat(1:(end-1),:),2)); % Vector of probabilities (areas)
            
            % Choose a region from which to sample
            aux=discreternd(cumsum(P),1);
            u{i}(l,:)=mod(ceil(aux./(2.^(0:(q+p-1))))+1,2);% Chosen region
            
            % Sample uniformly in chosen region
            U{i}(l,:)=rand(1,q+p).*exp(lP_mat(aux,:))+(1-exp(lP_mat(aux,:))).*(u{i}(l,:)==0);
            if sum(isnan(U{i}(l,:)))>0
                u{i}(l,:)
                U{i}(l,:)
            end
        end
    end
end