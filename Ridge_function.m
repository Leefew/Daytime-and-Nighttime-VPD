function [cof,r,pvl,k]=myRR_vif(Y,X,vif_limit,k_setp,k_limit)
[n,p]=size(X);
mx=mean(X);
my=mean(Y);
stdx=std(X);
stdy=std(Y);
MX=mx(ones(n,1),:);
STDX=stdx(ones(n,1),:);
Z=(X-MX)./STDX;
Y=(Y-my)./stdy;

kz=0;
while kz<=k_limit
    pseudo=sqrt(kz*(n-1))*eye(p);
    Zplus=[Z;pseudo];%;
    
    rc=corrcoef(Zplus);
    %                 vi=cond(rc);
    vi=max(diag(inv(rc)));
    if vi<vif_limit
        break
    end
    kz=kz+k_setp;
    
end

pseudo2=sqrt(kz*(n-1))*eye(p);
% Zplus=[Z;pseudo];%;ones(1,size(Z,2))
Zplus=[Z;pseudo2];%;
Zplus2=[Zplus,ones(size(Zplus,1),1)];
Yplus=[Y;zeros(p,1)];
% [b,bint,r,rint,stats]=regress(Yplus,Zplus2);
[b,~,~,~,stats]=regress(Yplus,Zplus2);
cof=b(1:end-1);
r=stats(1);
pvl=stats(3);
k=kz;
end

