function alist = comm_sync(Aset,M,m)

M = M(:); m = m(:);
[~,~,sims] = size(Aset);
nC = sum(m);
nS = numel(M);
urows = [0;cumsum(M.*m)];

alist = zeros(nC,sims);
cc = 0;
for i=1:nS
  for c = 1:m(i)
    cc = cc + 1;
    srow = urows(i)+(c-1)*M(i)+1;
    erow = urows(i)+c*M(i);
    alist(cc,:) = squeeze(mean(mean(Aset(srow:erow,srow:erow,:),1),2));
  end
end
