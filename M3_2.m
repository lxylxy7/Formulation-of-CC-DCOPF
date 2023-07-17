clear
%%


load('D:\CC_DCOPF程序\case118_cosample_LW_1000_2.mat');
load('D:\CC_DCOPF程序\case118_GMM_1000_15')

num=1000;
mpc=case118;
mpc2=case118WashingtonData;
mpc.branch(:,6)=mpc2.branch(:,6)*1.5;
% mpc.branch(:,6)=100000;
mpc.gen(:,9)=mpc.gen(:,9)*1.5;

bus_num=size(mpc.bus,1);
branch_num=size(mpc.branch,1);
gen_num=size(mpc.gen,1);
W_num=0;
P_max=mpc.gen(:,9);
P_min=mpc.gen(:,10);

branch_max=mpc.branch(:,6);
% branch_max(2:5)=100000;

PD=mu_D-W;

%%
b=1./mpc.branch(:,4);%电抗求倒数成电纳
Bbus=zeros(bus_num,bus_num);
f=mpc.branch(:,1);
t=mpc.branch(:,2);
nl = [(1:branch_num)'; (1:branch_num)'];     
Cft = sparse(nl, [f;t], [ones(branch_num, 1); -ones(branch_num, 1)],branch_num , bus_num);
Bf = sparse(nl, [f; t], [b; -b]); %branch_num x bus_num
Bbus = Cft' * Bf; %  bus_num x bus_num
power_gen=mpc.gen(:,1);%发电机对应节点
Ag=sparse(power_gen,1:gen_num,ones(gen_num,1),bus_num,gen_num); %火电机组对应节点关系

slack_bus = 1;   %平衡节点设置
noref   = (2:bus_num)';      %将1作为电压参考节点
noslack = find((1:bus_num)' ~= slack_bus);
PTDF= zeros(branch_num, bus_num);
PTDF(:, noslack) = full(Bf(:, noref) / Bbus(noslack, noref));
Beta=mpc.gen(:,9)/sum(mpc.gen(:,9));


PTDF=makePTDF(mpc.baseMVA,mpc.bus,mpc.branch,4);

%%
p=sdpvar(gen_num,1,'full');

%% 分段线性所需材料

%%
fmin=[mpc.gencost(:,5)*100+mpc.gencost(:,6)]'*p;

%%
% cons=[];
lv=[sum(PD)-sum(p)==0];
branch_cons=[];
for  br=1:branch_num
   
branch_cons=[branch_cons,PTDF(br,:)*(Ag*p-PD)<=branch_max(br)-PF_dn_gm(br)]; %
branch_cons=[branch_cons,PTDF(br,:)*(Ag*p-PD)>=-branch_max(br)+PF_up_gm(br)]; %
end

gen_cons=[];
for i=1:gen_num
    
gen_cons=[gen_cons,p(i)<=P_max(i)-G_up_gm(i)];
gen_cons=[gen_cons,p(i)>=P_min(i)+G_dn_gm(i)];%
end
% b_cons=[];
% b_cons=[b_cons,Beta>=0];
% b_cons=[b_cons,sum(Beta)==1];
cons=[lv,branch_cons,gen_cons];


%%

ops=sdpsettings('solver', 'gurobi','savesolveroutput',1);
tic
result=optimize(cons,fmin,ops);
solution_time=toc
Afmin=value(fmin) ;

app=value(p);
p=value(p);
Beta=value(Beta);
% UM_dn=value(UM_dn);
% % UM_up=value(UM_up);
% KK=value(KK);
% K=value(K);
% AA1=value(AA1);
% AA2=value(AA2);
% AA3=value(APFF);
% APPP2=value(APPP2);
% load('D:\程序\威布尔\Sample_W.mat')




% AAp=-Beta.*sum(WW-W,1);
% p=W+AAp;

DD=bus_p-bus_w;
AAp=Beta.*sum(DD-PD,1);
p=p+AAp;
k=0;

br=PTDF*(Ag*p-DD);

% Ab6_dn=size(find(br(6,:)<-240),2)/num;
for i=1:size(br,1)
  Ab_up(i)=size(find(br(i,:)>branch_max(i)),2)/num; 
  Ab_dn(i)=size(find(br(i,:)<-branch_max(i)-0.1),2)/num;
end


for i=1:size(p,1)
Ag_up(i)=size(find(p(i,:)>P_max(i)),2)/num;
Ag_dn(i)=size(find(p(i,:)<0),2)/num;

end
Abr=max([Ab_up,Ab_dn]);
Agg=max([Ag_up,Ag_dn]);
