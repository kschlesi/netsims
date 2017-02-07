function [theta_ens,A_ens,C_ens] = ksims_ens(sims,struct_gen,varargin)
%%% runs an ensemble of sims with ksims function %%%
% inputs: 'sims' - number of sims to run
%         'struct_gen' - function taking 0 args and returning struct matrix 
%         'varargin' - (kappa_,sigma_,ts,endtime) and other args for ksims
% outputs: 'A_ens' - the NxNxts synchronization matrix for all sims
%          'C_ens' - the NxN matrix for all sims
%          'theta_ens' - the Nxts matrix of thetas for all sims

% initialize output matrices
ts = varargin{3};
endtime = varargin{4};
[N,~] = size(struct_gen());
A_ens = zeros(numel(0:ts:endtime),N,N,sims);
C_ens = zeros(N,N,sims);
theta_ens = zeros(numel(0:ts:endtime),N,sims);

% loop over all simulations
for s=1:sims
    % create stochastic struct matrix
    C = struct_gen();
    
    % run simulation
    [theta,A] = ksims(C,varargin{:});
    
    % save C, A, and theta in return matrices
    A_ens(:,:,:,s) = A;
    C_ens(:,:,s) = C;
    theta_ens(:,:,s) = theta;
    
end

end