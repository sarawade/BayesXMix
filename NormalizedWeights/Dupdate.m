function [da,Da]=Dupdate(phi,Y,x,n,q,p,sigmaa,betaa,rhoa,mua,taua,wa,Jd,ka,JD)
% function [d,D]=Dupdate(phi,Y,x,n,q,p,sigmaa,betaa,rhoa,mua,taua,wa,Jd,ka,JD)
% Updates the index latent variables with slice sampling techniques
% INPUT:
% phi: Slice sampling constant
% Y: n x 1 continuous response variable
% x: n x q+p matrix of covariates 
% n: sample size
% q: number of discrete covariates 
% p: number of continuous covariates 
% betaa: q+p+1 x Ja matrix of current regression coefficients
% sigmaa: 1 x Ja vector of current regression variance parameters
% rhoa: cell array of length q, where cell h is a matrix of size Gh x Ja of
%   current discrete covariate parameters
% mua: p x Ja matrix of current mean parameters
% taua: p x 1 vector with the current state of the common component
%   precision
% wa: vector of length Ja of current mixture weights
% Jd: vector of number of components needed to sample the d
% ka: 1 x n vector with current k values
% JD: cell array of number of components needed to sample the D
%
% OUTPUT:
% da: 1 x n vector with updated state of index variables
% Da: cell array of length n with updated index values for latent variables

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ALLOCATION %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
da=zeros(1,n);
Da=cell(n,1);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAMPLE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:n
%%  Indices associated to the data (d)
    if Jd(i)==1 % Only one component available
        da(i)=1;
    else % More than one component is possible
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CALCULATE POSTERIOR PROBABILITIES %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        lpp=zeros(Jd(i),1); % Initialize log posterior probabilities
        for j=1:Jd(i)
            % Compute log likelihood
            lg1=0;
            lg2=0;
            if q>0 % Discrete covariate contribution
                for h=1:q
                    Gh=size(rhoa{h},1);
                    lg1=lg1+sum(log(rhoa{h}(:,j)).*(x(i,h)==(0:Gh-1)'));
                end
            end
            if p>0 % Continuous variable contribution
                lg2=-0.5*sum((x(i,q+1:end)'-mua(:,j)).^2.*taua);
            end
            % Response variable contribution
            lgy=-0.5*(Y(i)-[1,x(i,:)]*betaa(:,j))^2/sigmaa(j);
            lpp(j)=log(wa(j))+j*phi+lg1+lg2+lgy-log(sigmaa(j))*.5;
        end
        % Find normalization constant
        funky = @(const) sum(exp(const+lpp))-1;
        [const,~]=fzero(funky,100);
        pp=exp(const+lpp); % Posterior probabilities
        % Sample d
        da(i)=discreternd(cumsum(pp),1);
    end
%%  Indices associated with latent variables (D)
    if ka(i)>0 % Lantent variables are present
        for l=1:ka(i)
            if JD{i}(l)==1 % Only one component available
                Da{i}(1,l)=1;
            else % More than one component is possible
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             % CALCULATE POSTERIOR PROBABILITIES %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                lpp=zeros(JD{i}(l),1);% Initialize log posterior probabilities
                for j=1:JD{i}(l)
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
                    lpp(j)=log(wa(j))+j*phi+log(1-exp(lg1+lg2));
                end
                % Find normalization constant
                funky = @(const) sum(exp(const+lpp))-1;
                if ~isreal(funky(100))||~isfinite(funky(100))
                    lpp
                end
                [const,~]=fzero(funky,100);
                pp=exp(const+lpp);% Posterior probabilities
                % Sample D
                Da{i}(1,l)=discreternd(cumsum(pp),1);
            end
        end
    end
end
