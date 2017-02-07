%%% Simulator for Kuramoto Oscillator Example

sims = 5;
endtime = 10;
ts = 0.1;
period = 2*pi;

% define parameters
kappa_ = 0.2;   % constant coupling strength
sigma_ = 1;     % standard deviation of intrinsic frequencies (norm dist)
N = 100;        % number of nodes
M = [20;10];    % unique community sizes present
m = [2;4];      % number of communities of each size
pin = [0.9,0.6];% coupling probability for in-community nodes
pout = 0.01;    % coupling probability for out-of-community nodes
pbase = pout;   % coupling probability for all singletons

%%
% run sims 
[theta_ens,A_ens,C_ens] = ksims(sims,N,M,m,pin,pout,pbase,kappa_,sigma_,...
                                            ts,endtime);

% figure: underlying connectivity matrix
figure; bcolor(C); 

% figure: position v. time (example)
mtheta_ens = mod(theta_ens,period);
figure; bcolor(mtheta_ens(:,:,end)'./period); colorbar;

% figure: final synchronization (example)
figure; bcolor(squeeze(A_ens(end,:,:,end))); colorbar;

% figure: final synchronization (mean)
figure; bcolor(squeeze(mean(A_ens(end,:,:,:),4))); colorbar;

%%
% do a sweep over pins
pinlist = 0:0.2:1;
A_fin = zeros(N,N,numel(pinlist),numel(pinlist));
C_sample = zeros(N,N,numel(pinlist),numel(pinlist));
f=0;
for p1 = pinlist
    for p2 = pinlist
        f = f+1;
        [~,A_ens,C_ens] = ksims(sims,N,M,m,[p1;p2],pout,pbase,kappa_,sigma_,...
                                            ts,endtime);
        A_fin(:,:,pinlist==p1,pinlist==p2) = squeeze(mean(A_ens(end,:,:,:),4));
        C_sample(:,:,pinlist==p1,pinlist==p2) = C_ens(:,:,end) ;
        figure(1); hold on;
        subplot(numel(pinlist),numel(pinlist),f);
        bcolor(A_fin(:,:,pinlist==p1,pinlist==p2)); 
        title(['p1 = ' num2str(p1) ', p2 = ' num2str(p2)]);
        hold off;
        figure(2); hold on;
        subplot(numel(pinlist),numel(pinlist),f);
        bcolor(C_sample(:,:,pinlist==p1,pinlist==p2));
        title(['p1 = ' num2str(p1) ', p2 = ' num2str(p2)]);
        hold off;
    end
end

%%
% do a sweep over community sizes
pinlist = 0.2:0.2:1;
Mlist = [5,10,20,30,50];
siglist = [0.5,1,5,10,50];
As = zeros(numel(pinlist),numel(siglist),numel(Mlist));
Ad = zeros(numel(pinlist),numel(siglist),numel(Mlist));
Cs = zeros(N,N,numel(pinlist),numel(Mlist));
Aend = zeros(N,N,numel(pinlist),numel(siglist),numel(Mlist));
for M = Mlist
  for pin = pinlist
    for sigma_ = siglist
        m = floor(N/M);
        [~,A_ens,C_ens] = ksims(sims,N,M,m,pin,pout,pbase,...
                                kappa_,sigma_,ts,endtime);
        Aend(:,:,pinlist==pin,siglist==sigma_,Mlist==M) = squeeze(mean(A_ens(end,:,:,:),4));
        % compute mean total synchronization in each community (assumes all same size)
        As(pinlist==pin,siglist==sigma_,Mlist==M) = ...
            mean(mean(comm_sync(squeeze(A_ens(end,:,:,:)),M,m)));
        % compute sumtotal synchronization between non-community pairs
        Ad(pinlist==pin,siglist==sigma_,Mlist==M) = ...
            sum(sum(squeeze(mean(A_ens(end,:,:,:),4)))) - m*As(pinlist==pin,siglist==sigma_,Mlist==M);
    end
    Cs(:,:,pinlist==pin,Mlist==M) = C_ens(:,:,end);
  end
end

%%
% plot things from previous run
f2 = 0;
for M = Mlist
  f = 0;
  for pin = pinlist
    for sigma_ = siglist
        f = f+1;
        figure(find(Mlist==M)); hold on;
        subplot(numel(pinlist),numel(siglist),f);
        bcolor(Aend(:,:,pinlist==pin,siglist==sigma_,Mlist==M));
        title(['M = ' num2str(M) ', pin = ' num2str(pin) ', sig = ' num2str(sigma_)]); 
        %suptitle(['M = ' num2str(M)]); 
        hold off;
    end
    f2 = f2+1;
    figure(2*numel(Mlist)+1); hold on;
    subplot(numel(pinlist),numel(Mlist),f2);
    bcolor(Cs(:,:,pinlist==pin,Mlist==M)); 
    title(['M = ' num2str(M) ', pin = ' num2str(pin)]); hold off;
  end
  % plot average sync per community, for each M, pin v sig
  figure(numel(Mlist)+find(Mlist==M)); hold on;
  subplot(1,2,1); bcolor(As(:,:,Mlist==M)');
  xlabel('in-comm connection prob'); ylabel('natural freq variance');
  set(gca,'XTick',(0:1:numel(pinlist)-1)+1.5,'YTick',(0:1:numel(siglist)-1)+1.5);
  set(gca,'XTickLabel',pinlist,'YTickLabel',siglist);
  title(['mean comm sync, M = ' num2str(M)]); hold off;
  subplot(1,2,2); bcolor(Ad(:,:,Mlist==M)');
  xlabel('in-comm connection prob'); ylabel('natural freq variance');
  set(gca,'XTick',(0:1:numel(pinlist)-1)+1.5,'YTick',(0:1:numel(siglist)-1)+1.5);
  set(gca,'XTickLabel',pinlist,'YTickLabel',siglist);
  title(['mean noncomm sync, M = ' num2str(M)]); hold off;
end
%%
%%%%%%%% STATIC community detection %%%%%%%%%
p = 100; gamma = 1; omega = 0;
Mlist = Mlist(1:3);
siglist = siglist(1:3);
C = zeros(p,N,numel(pinlist),numel(siglist),numel(Mlist));
Cnew = zeros(p,N,numel(pinlist),numel(siglist),numel(Mlist));
Q = zeros(p,numel(pinlist),numel(siglist),numel(Mlist));
Qnew = zeros(p,numel(pinlist),numel(siglist),numel(Mlist));
for M=Mlist
  for pin = pinlist
     for sigma_ = siglist
       % static community detection
       [C(:,:,pinlist==pin,siglist==sigma_,Mlist==M),Q(:,pinlist==pin,siglist==sigma_,Mlist==M)] = ...
           genlouvainREPs(Aend(:,:,pinlist==pin,siglist==sigma_,Mlist==M),p,gamma,omega);
       % community consensus
       [Cnew(:,:,pinlist==pin,siglist==sigma_,Mlist==M),Qnew(:,pinlist==pin,siglist==sigma_,Mlist==M)] = ...
           consensus_comm_GL2(C(:,:,pinlist==pin,siglist==sigma_,Mlist==M));  % uses "genlouvainREPs"
     end
  end
end


%%
% community detection results

for M = Mlist
  f = 0;
  f1 = 0;
  for pin = pinlist
    for sigma_ = siglist
      % plot original communities
      f = f+1;
      figure(find(Mlist==M)); hold on;
      subplot(numel(pinlist),numel(siglist),f);
      bcolor(C(:,:,pinlist==pin,siglist==sigma_,Mlist==M)');
      title(['M = ' num2str(M) ', pin = ' num2str(pin) ', sig = ' num2str(sigma_)]);
      hold off;
      % plot CC communities
      f1 = f1+1;
      figure(find(Mlist==M)+numel(Mlist)); hold on;
      subplot(numel(pinlist),numel(siglist),f1);
      bcolor(Cnew(:,:,pinlist==pin,siglist==sigma_,Mlist==M)');
      title(['M = ' num2str(M) ', pin = ' num2str(pin) ', sig = ' num2str(sigma_)]);
      hold off;
    end
  end
  % make a sig v pin plot of Q's
  figure(2*numel(Mlist)+find(Mlist==M)); hold on;
  bcolor(squeeze(mean(Q(:,:,:,Mlist==M),1))');
  xlabel('in-comm connection prob'); ylabel('natural freq variance');
  set(gca,'XTick',(0:1:numel(pinlist)-1)+1.5,'YTick',(0:1:numel(siglist)-1)+1.5);
  set(gca,'XTickLabel',pinlist,'YTickLabel',siglist);
  title(['mean Q, M = ' num2str(M)]); colorbar; hold off;
end

%%
%%%%%%%% plots %%%%%%%%%

% from single realisation: phase of each oscillator over time
figure
plot(tspan,theta(:,1))
hold on
hold all
for i=2:N
    plot(tspan,theta(:,i))
end

% single realisation: final phase of each oscillator
mtheta = mod(theta,period);
figure
plot(1:N,sort(mtheta(end,:))./period)

% synchronization of osc 1 with all others over time (avg of sims realzns)
figure
plot(tspan,A(:,:,1))

