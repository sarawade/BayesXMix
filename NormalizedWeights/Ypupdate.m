function [betaa,sigmaa]=Ypupdate(Y,X,da,Ja,beta0,iC,alpha1,alpha2,C)
% function [betaa,sigmaa]=Ypupdate(Y,X,da,Ja,beta0,C,alpha1,alpha2)
% Updates the regression parameters (beta, sigma) with Normal kernel and
% multivariate Normal-InverseGamma base measure
% INPUT:
% Y: n x 1 continuous response variable
% X: n x q+p+1 augmented covariate matrix
% da: 1 x n vector with current state of index variables
% Ja: current number of components
% beta0,iC,alpha1,alpha2: hyperparameters for the base measure
% C=(iC)^-1
%
% OUTPUT:
% betaa: q+p+1 x Ja matrix of updated regression coefficients
% sigmaa: 1 x Ja vector of updated regression variance parameters
  
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ALLOCATION %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
betaa=zeros(length(beta0),Ja);
sigmaa= zeros(1,Ja);
 	
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAMPLE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 for j=1:Ja
 	n_j=sum(da==j); % number of observations associated to component j
    if n_j>0 % There are observations associated to component j
        X_j=X(da==j,:); % X matrix of covariates associated to component j
        Y_j=Y(da==j); % Response variables associated to component j
        % Posterior parameters for beta
        C_hat=C+X_j'*X_j; % posterior precision matrix
        % The matrix is positive semi-definite, so we use eigen value decomposition to compute inverse
        [V,lambda]=eig(C_hat);
        iC_hat=V*diag(1./diag(lambda))*V';
        beta0_hat=iC_hat*(C*beta0+X_j'*Y_j); % posterior mean
        % Posterior parameters for sigma
        alpha1_hat=alpha1+n_j/2;
        W_j=eye(n_j)-X_j*iC_hat*X_j';
        alpha2_hat=alpha2+1/2*(Y_j-X_j*beta0)'*W_j*(Y_j-X_j*beta0);
    else % Otherwise, sample from the prior
        iC_hat=iC;
        beta0_hat=beta0;
        alpha1_hat=alpha1;
        alpha2_hat=alpha2;
    end
    % Sample
    sigmaa(j)=1/gamrnd(alpha1_hat,1/ alpha2_hat);
    betaa(:,j)=chol(iC_hat,'lower')*randn(length(beta0),1)*sigmaa(j)^.5+beta0_hat;
end

