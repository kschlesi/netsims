function [theta,A] = ksims(C,kappa_,sigma_,ts,endtime)
%%%%%%% run simulation %%%%%%%%

% read in structure matrix
try
    assert(ismatrix(C));
    [N,N1] = size(C);
    assert(N1==N);
catch
    error('Input structure matrix must be 2D and square.');
end

% initialize simulation time
tspan = 0:ts:endtime;

% create initial condition vectors
w = randn(N,1).*sigma_;                % frequencies (normally distributed)

% initialize state variables
theta0 = ones(N,1);
A = zeros(numel(tspan),N,N);

% solve discretised equation & keep results in theta
theta = zeros(numel(tspan),numel(theta0));
theta(1,:) = theta0;
for t=2:numel(tspan)
  
  % evolve oscillators forward in time
  sinji = sin(repmat(theta(t-1,:),N,1)-repmat(theta(t-1,:)',1,N));
  theta(t,:) = theta(t-1,:)' + ts.*w + kappa_.*sum(C.*sinji,2);      
   
  % compute synchronization matrix
  cosij = cos(repmat(theta(t,:)',1,N)-repmat(theta(t,:),N,1));
  %A(t,:,:) = A(t,:,:) + shiftdim(abs(cosij)./sims,-1);
  A(t,:,:) = shiftdim(abs(cosij),-1);
  
end  % end solve over tspan

end