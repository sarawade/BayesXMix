function [X]=tgamrnd(type,B,a,b,X,burnin)
% function [X]=tgamrnd(type,B,a,b,X,burnin)
% Samples from a truncated Gamma(a,b) distribution.
% Using Gibbs sampler and latent variable technique similar to that of Damien, Walker
% (2011)Sampling Truncated Normal, Beta, and Gamma Densities Journal of 
% Computational and Graphical Statistics, Vol. 10, No. 2, pp. 206-215
% 
% INPUT:
%   type: may be
%       'left' if upper bound, so sample form left tail
%       'right' if lower bound, so sample from right tail
%       'between' if sampling from a bounded interval
%       B: upper bound if type='left'
%          lower bound if type='right'
%          [lbound,rbound] if type='between'
%   (a,b): parameters for the Gamma with mean a/b
%   X: initial value for MCMC simulation 
%   burnin: number of iterations for the burn in period
%
% OUTPUT: 
%   X= one realization of the truncated Gamma.

if strcmp(type,'right') % Sampling from right tail
    for i=1:burnin+1
        B_v=0;
        if a~=1
            B_v=exp(log(rand)*(1/(a-1))+log(X));
        end
        B_v=max(B,B_v);
        X=exprnd(1/b)+B_v;
    end
else
    if strcmp(type,'between') % Sampling from bounded interval
        B_1=B(1);
        B_2=B(2);
    else % Sampling from left tail
        B_1=0;
        B_2=B;
    end
    for i=1:burnin+1
        B_v=0;
        if a~=1
            B_v=exp(log(rand)*(1/(a-1))+log(X));
        end
        B_v=max(B_1,B_v);
        I=1-exp(-b*(B_2-B_v));
        if I==0 % Precision issues; approximation required
            disp('I am approximating truncated gamma!')
            X=B_v+(B_2-B_v)*rand;
        else
            X=B_v-1/b*log(1-rand*I);
        end
    end
end
