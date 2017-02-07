function [C,Q,C_hist,npasses] = genlouvainREPs(A,p,gamma,omega,max_passes)
% this code does p community detection optimizations on adjacency matrix A
% using generalized Louvain algorithm from netwiki package.
% genlouvain returns S, an nx1 vector of assignments, and a 1x1 Q-value.

% inputs: A (original NxN (potentially sparse) adj matrix for static)
%           OR (Tx1 cell array of T NxN matrices for categorical multislice)
%         gamma (value of gamma for use in community optimization penalty)
%         omega (value of omega for multislice; omit or set = 0 for static)
%         p (number of genlouvain optimizations to perform)
%         max_passes (max # of passes whose history to save; default = 7)
%
% outputs: C (pxN matrix of community assignments)
%            OR (pxNxT tensor of multislice community assignments)
%          Q (px1 matrix of associated quality values)
%          C_hist (pxNxTx(max_npasses) tensor of iterative C values)

if nargin<5 % set default max_passes
    max_passes = 7;
end

if nargin<4 % set default omega (static)
    omega = 0;
end

if ~iscell(A) % A is an NxN matrix
    A = {A};
end


% multislice community detection, A is Tx1 array of NxN matrices
    N=length(A{1});
    T=length(A);
    B=spalloc(N*T,N*T,(N+T)*N*T);
    twomu=0;
    for s=1:T
        k=sum(A{s});
        twom=sum(k);
        twomu=twomu+twom;
        indx=(1:N)+(s-1)*N;
        B(indx,indx)=A{s}-gamma*(k'*k)/twom;
    end
    twomu=twomu+T*omega*N*(T-1);
    all2all = N*[(-T+1):-1,1:(T-1)];
    B = B + omega*spdiags(ones(N*T,2*T-2),all2all,N*T,N*T);
    C = zeros(p,N,T);
    Q = zeros(p,1);
    C_hist = zeros(p,N,T,max_passes);
    npasses = zeros(p,1);
    for i=1:p  % performing p optimizations
        disp(i);
        %[S,Q1,S_hist] = genlouvainmex(B);
        [S,Q1,S_hist] = genlouvain_hist(B);
        Q(i) = Q1/twomu;
        C(i,:,:) = reshape(S,N,T);
        npasses(i) = size(S_hist,2);
        S_hist = reshape(S_hist,[N,T,size(S_hist,2)]);
        C_hist(i,:,:,1:npasses(i)) = shiftdim(S_hist,-1);
    end
    C_hist = C_hist(:,:,:,1:max(npasses));
end