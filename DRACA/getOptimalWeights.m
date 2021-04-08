function [p,Ws]=getOptimalWeights(Hs,alpha)
%% seek p and compute Ws, see Algorithm 1 for more detail
% Hs --- a R x 1 vector, the r-th element is the trace of |R-GSGt|^2 on the r-th kernel
% alpha --- the tradeoff parameter on W
K=length(Hs);
[sorted,index]=sort(Hs,'ascend');
newHs=Hs(index);
p=K;
bfind=1;
gamma=0;
while p>0 && bfind
     gamma=(sum(newHs(1:p))+2*alpha)/p;
     if (gamma-newHs(p))>0
         bfind=0;
     else
         p=p-1;
     end
end
newWs=zeros(K,1);
for ii=1:p
      newWs(ii)=(gamma-newHs(ii))/(2*alpha);
end
Ws=zeros(K,1);
Ws(index)=newWs;
