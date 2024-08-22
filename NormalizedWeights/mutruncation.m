function [T,B]=mutruncation(j,h,x,n,q,tauh,ka,da,U,u)
% function [T,B]=mutruncation(j,h,x,n,q,taua,ka,da,U,u)
% Finds the truncation for mu(h,j) 
% INPUT:
% j: current component
% h: current continuous covariate coordinate
% x: n x 1 vector of continuous covariate h 
% n: sample size
% q: number of discrete covariates 
% tauh: current precition for continuous covariate h
% ka: 1 x n vector with current k values
% da: 1 x n vector with current state of index variables
% (U,u): cell arrays of length n with current uniform latent variables and
%   associated indicators
%
% OUTPUT:
% T: Indicator of posterior support including the tails
% B: Vector containing boundaries of intervals of posterior support

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B=[]; 
T=1; 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIND TRUNCATION INTERVALS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:n
    if ka(i)>0 % There is a truncation region associated with x(i)
        % Find the bounds for the regions associated with x(i) 
        s0=0; 
        s1=0;
        g0_ind=0; % Indicator of u_ilh = 0 for some l
        g1_ind=0; % Indicator of u_ilh = 1 for some l 
        for l=1:ka(i)
            if da{i}(l)==j
                % Calculate size of interval around x(i)
                s=sqrt((-2*log(U{i}(l,h+q)))/tauh);
                if u{i}(l,h+q)==0 % Support excludes interval around x(i)
                    % Calculate largest interval around x(i)
                    s0=max(s0,s);
                    g0_ind=1;
                else % Support excludes the complement of interval around x(i)
                    % Calculate smallest interval around x(i)
                    s1=min(s1,s)*g1_ind+s*(g1_ind==0);
                    g1_ind=1;
                end
            end
        end
        if g1_ind==1 || g0_ind==1 % x(i) defines a region for mu_j
            % There are 3 possible shapes for the region
            Ti=(g1_ind==0); % Indicator of this region including the tails
            if g0_ind*g1_ind==1 % Region is two bounded intervals
                Bi=[x(i)-s1,x(i)-s0,x(i)+s0,x(i)+s1];
            elseif g1_ind==1 % Region is one bounded interval
                Bi=[x(i)-s1,x(i)+s1];
            else % Region is two tails
                Bi=[x(i)-s0,x(i)+s0];
            end
        
            %% Take Intersection with previous regions in B 
            % Current value of B is defined by the intersection of regions associated with x(i'), i'<i
            if isempty(B) % No previous regions defined
                B=Bi;
                T=Ti;
            else % Intersection must be found
                % Find relative position of the Bi elements with respect to
                % the B elements
                Bi_pos=sum(kron(ones(size(Bi,2),1),B)<kron(ones(1,size(B,2)),Bi'),2);
                % Eliminate unnecesary boundary points from B
                if g0_ind*g1_ind==1 % Eliminate boundary points outside two bounded intervals
                    B_new=B([(Bi_pos(1)+1):Bi_pos(2),(Bi_pos(3)+1):Bi_pos(4)]); 
                elseif g1_ind==1 % Eliminate boundary points outside one bounded interval
                    B_new=B((Bi_pos(1)+1):Bi_pos(2));
                else % Eliminate boundary points outside of two tails
                    B_new=B([1:(Bi_pos(1)),(Bi_pos(2)+1):size(B,2)]);
                end
                % Indicator of the elements of Bi which need to be added to
                % B 
                toadd=((mod(Bi_pos,2)==1).*(T==0)+(mod(Bi_pos+1,2)==1).*(T==1)); 
                    %If B includes the tails, the elements of Bi with odd
                        % relative position are added; 
                    %otherwise, the elements with even relative position are added.
                B=sort([B_new, Bi(toadd==1)]); % Add elements and reorder B
                T=min(T,Ti); % Reasess if B includes the tails
            end 
        end
    end
end