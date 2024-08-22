function [ka,Da]=kupdate2(pk,x,n,q,p,Ja,wa,rhoa,taua, mua, ka, Da)
% function [ka,Da]=kupdate(pk,x,n,q,p,Ja,wa,rhoa,taua, mua, ka, Da)
% Updates the latent k variables
% INPUT:
% pk: Reversible jump constant
% x: n x q+p matrix of covariates 
% n: sample size
% q: number of discrete covariates 
% p: number of continuous covariates 
% wa: vector of length Ja of current mixture weights
% rhoa: cell array of length q, where cell h is a matrix of size Gh x Ja of
%   current discrete covariate parameters
% taua: p x 1 vector with the current state of the common component
%   precision
% mua: p x Ja matrix of current mean parameters
% ka: 1 x n vector with current k values
% Da: cell array of length n with current index values for latent variables
%
% OUTPUT:
% ka: 1 x n vector with updated k values
% Da: cell array of length n with updated index values for latent variables

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REVERSIBLE JUMP UPDATE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:n
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CALCULATE POSTERIOR PROBABILITIES %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lpp=zeros(Ja,1);% Initialize log posterior probabilities
    for j=1:Ja
        % Compute log likelihood
        lg1=0;
        lg2=0;
        if q>0 % Discrete covariate contribution
            for h=1:q
                Gh=size(rhoa{h},1);
                lg1=lg1+sum(log(rhoa{h}(:,j)).*(x(i,h)==(0:Gh-1)'));
            end
        end
        if p>0 % Continuous covariate contribution
            lg2=-0.5*sum((x(i,q+1:end)'-mua(:,j)).^2.*taua);
        end
            lpp(j)=log(wa(j))+log(1-exp(lg1+lg2));
    end
    % Find normalization constant
    funky = @(const) sum(exp(const+lpp))-1;
    [const,~]=fzero(funky,100);
    pp=exp(const+lpp);% Posterior probabilities
    % Sample D
    %% Moving from k up to k+1
    if rand<pk || ka(i)==0 % Propose the move
        % Sample new index
        Dnew=discreternd(cumsum(pp),1);
        % Calculate acceptance probability
        paccept=(1-pk)*exp(-const)/(pk);
        if rand<paccept % Accept the move
            ka(i)=ka(i)+1; % New k
            Da{i}(ka(i))=Dnew; % Add new index
        end
    else
        %% Moving from k down to k-1
        % Calculate acceptance probability
        paccept=pk/((1-pk))*exp(const);
        if rand<paccept % Accept the move
            ka(i)=ka(i)-1; % New k
            Da{i}= Da{i}(1:ka(i)); % Remove last index
        end
    end
end