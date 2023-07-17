clear
load('D:\CC_DCOPF程序\case118_cosample_LW_1000_2.mat');
U_bus_p=bus_p-mu_D;
U_bus_w=bus_w-W;
U_bus_p=U_bus_p-U_bus_w;
num=1000;

eeeee=0.15;
PD=mu_D-W;
mpc=case118;
mpc2=case118WashingtonData;
mpc.branch(:,6)=mpc2.branch(:,6)*1.5;
% mpc.branch(:,6)=100000;
mpc.gen(:,9)=mpc.gen(:,9)*1.5;
% kk=rundcopf(mpc);
bus_num=size(mpc.bus,1);
branch_num=size(mpc.branch,1);
gen_num=size(mpc.gen,1);
W_num=3;
P_max=mpc.gen(:,9);
P_min=mpc.gen(:,10);


D_bus=find(PD~=0);
D_num=size(D_bus,1);

branch_max=mpc.branch(:,6);
% branch_max(2:5)=100000;

power_gen=mpc.gen(:,1);%发电机对应节点
Ag=sparse(power_gen,1:gen_num,ones(gen_num,1),bus_num,gen_num); 
Ag_d=sparse(D_bus,1:D_num,ones(D_num,1),bus_num,D_num); %火电机组对应节点关系
PTDF=makePTDF(mpc.baseMVA,mpc.bus,mpc.branch,4);
Beta=mpc.gen(:,9)/sum(mpc.gen(:,9));
%% mc
for i=1:bus_num
   K(:,i)=-PTDF*Ag*Beta+PTDF(:,i); 
end
PF_UM_A=K*U_bus_p;% 正负号需要甄别
G_UM_A=Beta*sum(U_bus_p);
PF_UM_As=sort(PF_UM_A,2);
G_UM_As=sort(G_UM_A,2);
PF_dn_A=-PF_UM_As(:,num*eeeee);
PF_up_A=PF_UM_As(:,num*(1-eeeee));
G_dn_A=-G_UM_As(:,num*eeeee);
G_up_A=G_UM_As(:,num*(1-eeeee));
%% gm

U_bus_p2=U_bus_p;
U_bus_p2(find(U_bus_p(:,1)==0),:)=[];
% load_num=size(U_bus_p2,1);
numComponents=20;
% options = statset('MaxIter',500,'Display','final');
% gm = fitgmdist(U_bus_p2',numComponents,'Options',options);
[~,gm,llh] = mixGaussEm(U_bus_p2,numComponents);
% [~,gm,llh] = mixGaussVb(U_bus_p2,10);
% numComponents=size(gm.m,2);
% weight=exp(gm.logW);
weight=gm.w;

% mu1=gm.mu';
% sigma1=gm.Sigma;
mu=Ag_d*gm.mu;

for i=1:numComponents
sigma(:,:,i)=Ag_d*gm.Sigma(:,:,i)*Ag_d';
end
% weight=gm.w;

syms x; 
  h4=0.00709;
  h3=-0.02918;
  h2=-0.06211;
  h1=0.424;
  h0=-0.0005;
for k=1:gen_num
for i=1:size(mu,2)
%   ym1=(x-Beta(k)*ones(1,3)*mu1(:,i))/norm([Beta(k)*ones(1,3)]*blkdiag(sigma1(:,:,i)^(1/2)),2);
  ym=(x-Beta(k)*ones(1,bus_num)*mu(:,i))/norm([Beta(k)*ones(1,bus_num)]*blkdiag(sigma(:,:,i)^(1/2)),2);
%   Beta(k)*ones(1,3)*mu1(:,i)
%   Beta(k)*ones(1,gen_num)*mu(:,i)
%   norm([Beta(k)*ones(1,3)]*blkdiag(sigma1(:,:,i)^(1/2)),2)
%   norm([Beta(k)*ones(1,gen_num)]*blkdiag(sigma(:,:,i)^(1/2)),2)
  fup(i)=0.5+h4*ym.^4+h3*ym.^3+h2*ym.^2+h1*ym+h0;
  fdn(i)=0.5-h4*ym.^4+h3*ym.^3-h2*ym.^2+h1*ym-h0;  
end
Fup=weight*fup';
Fdn=weight*fdn';
% f=x^3-2*x+1;
Sup=solve(Fup==(1-eeeee),x,'Real',true);
Sdn=solve(Fdn==eeeee,x,'Real',true);
Sup=double(Sup);
Sdn=-double(Sdn);
for i=1:2
    if Sup(i)>=0
G_up_gm(k,:)=Sup(i);
    end
    if Sdn(i)>=0
G_dn_gm(k,:)=Sdn(i);
    end
end
end
clear fup fdn Fup Fdn Sup Sdn
for k=1:branch_num
for i=1:size(mu,2)

  ym=(x-K(k,:)*mu(:,i))/norm(K(k,:)*blkdiag(sigma(:,:,i)^(1/2)),2);
%   if K(k,:)*mu(:,i)>=0
  fup(i)=0.5+h4*ym.^4+h3*ym.^3+h2*ym.^2+h1*ym+h0;
  fdn(i)=0.5-h4*ym.^4+h3*ym.^3-h2*ym.^2+h1*ym-h0;   
%   else
%   fdn(i)=0.5+h4*ym.^4+h3*ym.^3+h2*ym.^2+h1*ym+h0;
%   fup(i)=0.5-h4*ym.^4+h3*ym.^3-h2*ym.^2+h1*ym-h0;      
%   end
end
Fup=weight*fup';
Fdn=weight*fdn';
% f=x^3-2*x+1;
Sup=solve(Fup==(1-eeeee),x,'Real',true);
Sdn=solve(Fdn==eeeee,x,'Real',true);
Sup=double(Sup);
Sdn=-double(Sdn);
for i=1:2
    if Sup(i)>=0
PF_up_gm(k,:)=Sup(i);
    end
    if Sdn(i)>=0
PF_dn_gm(k,:)=Sdn(i);
    end
end
end
% PF_up_gm=PF_up_A;
% PF_dn_gm=PF_dn_A;

PFUP=[PF_up_A,PF_up_gm];
PFDN=[PF_dn_A,PF_dn_gm];
save('D:\CC_DCOPF程序\case118_GMM_1000_15','PF_up_gm','PF_dn_gm','G_up_gm','G_dn_gm','PF_up_A','PF_dn_A','G_up_A','G_dn_A')