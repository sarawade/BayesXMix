function [ka,Da]=kupdate(pk,x,n,q,p,wa,rhoa,taua, mua, ka, Da)
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
cw=cumsum(wa); % Cumulative weights   

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REVERSIBLE JUMP UPDATE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:n
    %% Moving from k up to k+1
    if rand<pk || ka(i)==0 % Propose the move
        Dnew=discreternd(cw,1); % Sample new index
        % Calculate acceptance probability
        lg1new=0;
        lg2new=0;
        if q>0 % Discrete variable contribution
            for h=1:q
                Gh=size(rhoa{h},1);
                lg1new=lg1new+sum(log(rhoa{h}(:,Dnew)).*(x(i,h)==(0:Gh-1)'));
            end
        end
        if p>0 % Continuous variable contribution
            lg2new=-0.5*sum((x(i,q+1:end)'-mua(:,Dnew)).^2.*taua);
        end
        gnew=exp(lg1new+lg2new);
        paccept=(1-pk)*(1-gnew)*sum(wa)/(pk);
        if rand<paccept % Accept the move
            ka(i)=ka(i)+1; % New k
            Da{i}(ka(i))=Dnew; % Add new index
        end
    else
        %% Moving from k down to k-1
        % Calculate acceptance probability
        lg1k=0;
        lg2k=0;
        if q>0 % Discrete variable contribution
            for h=1:q
                Gh=size(rhoa{h},1);
                lg1k=lg1k+sum(log(rhoa{h}(:,Da{i}(ka(i)))).*(x(i,h)==(0:Gh-1)'));
            end
        end
        if p>0 % Continuous variable contribution
            lg2k=-0.5*sum((x(i,q+1:end)'-mua(:,Da{i}(ka(i)))).^2.*taua);
        end
        gk=exp(lg1k+lg2k);
        paccept=pk/((1-pk)*(1-gk)*sum(wa));
        if rand<paccept % Accept the move
            ka(i)=ka(i)-1; % New k
            Da{i}= Da{i}(1:ka(i)); % Remove last index
        end
    end
end