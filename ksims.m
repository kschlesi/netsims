function [theta_ens,A_ens,C_ens] = ksims(sims,N,M,m,pin,pout,pbase,kappa_,sigma_,...
                                            ts,endtime)
%%%%%%% run simulation %%%%%%%%
tspan = 0:ts:endtime;
A_ens = zeros(numel(tspan),N,N,sims);
C_ens = zeros(N,N,sims);
theta_ens = zeros(numel(tspan),N,sims);
for s=1:sims
    disp([num2str(s) ' of ' num2str(sims)]);

% create initial condition vectors
w = randn(N,1).*sigma_;                    % frequencies (normally distributed)
C = modcoupler(N,M,m,pbase,pin,pout);      % NxN binary coupling matrix
% i=0;
% while any(sum(C(1:M'*m,1:M'*m))==0)
%     disp('ensuring no disconnected in-community nodes');
%     disp(i); i=i+1;
%     C = modcoupler(N,M,m,pbase,pin,pout);  % NxN binary coupling matrix
% end

% initialize state variables
theta0 = ones(N,1);
A = zeros(numel(tspan),N,N);

% solve discretised equation & keep results in theta
theta = zeros(numel(tspan),numel(theta0));
for t=1:numel(tspan)
    
  if tspan(t)==0
      % initial condition
      theta(t,:) = theta0;
  else
      sinji = sin(repmat(theta(t-1,:),N,1)-repmat(theta(t-1,:)',1,N));
      theta(t,:) = theta(t-1,:)' + ts.*w + kappa_.*sum(C.*sinji,2);      
  end
   
  % compute synchronization matrix
  cosij = cos(repmat(theta(t,:)',1,N)-repmat(theta(t,:),N,1));
  %A(t,:,:) = A(t,:,:) + shiftdim(abs(cosij)./sims,-1);
  A_ens(t,:,:,s) = shiftdim(abs(cosij),-1);
  
end  % end solve over tspan

% save network struture and function from this instance
C_ens(:,:,s) = C;
theta_ens(:,:,s) = theta;

end  % end loop over simulations

end