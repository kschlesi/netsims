function [S2,Q2,X_new3,qpc] = consensus_comm_GL2(C,gamma,omega)

% FUNCTION: To identify consistent communities in an ensemble
% 
% OUTPUT:
% S2 is a pxn matrix of new community assignments
% Q2 is the associated modularity value
% X_new3 is the thresholded nodal association matrix
% qpc is the quality of the consensus (the lower the value the better)
%
% INPUT:
% C is a pxn OR pxnxt matrix of community assignments where p is the number of optimizations
% that have been performed. n is the number of nodes in the system. C gives the real partitions. 

% find size of network & number of optimizations
[p,n,t] = size(C);
reshape(C,p,n*t); % unfold into 2 dimensions

% compute thresholded module allegiance matrix with random permutations of C
X_new3 = mod_allegiance(C,1).*p;

% recompute optimal partition on this new matrix of kept community
% association assignments (output: S2, the new pxnxt assignment matrix)
X_new = cell(t,1);
for T=1:t
    X_new{T} = X_new3(n*(T-1)+1:n*T,n*(T-1)+1:n*T);
end
[S2,Q2] = genlouvainREPs_f(X_new,p,gamma,omega);

% define the quality of the consensus
qpc = sum(sum(abs(diff(S2))));