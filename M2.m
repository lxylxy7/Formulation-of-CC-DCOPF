clear

load('nearoptimalsolutions.mat')
load('case118_uncertaintydata_1000.mat');

%% near optimal solutions loading
DD=bus_p-bus_w;
AAp=Beta.*sum(DD-PD,1);
p=app+AAp;
br=PTDF*(Ag*p-DD);

PP=p;
BR=br;

BETA=Beta;
clear p br Beta



%% power system definition
num=1000; % scenario number
eeee=0.05;% predefined violation probability

bus_p=bus_p(:,1:num);
bus_w=bus_w(:,1:num);
% 
mpc=case118;
mpc2=case118WashingtonData;
mpc.branch(:,6)=mpc2.branch(:,6)*1.5;

mpc.gen(:,9)=mpc.gen(:,9)*1.5;


bus_num=size(mpc.bus,1);
branch_num=size(mpc.branch,1);
gen_num=size(mpc.gen,1);

P_max=mpc.gen(:,9);
P_min=mpc.gen(:,10);

branch_max=mpc.branch(:,6);

PD=mu_D-W; % forecast power value of buses

U_bus=bus_p-bus_w-PD;% uncertainty of buses

b=1./mpc.branch(:,4);
Bbus=zeros(bus_num,bus_num);
f=mpc.branch(:,1);
t=mpc.branch(:,2);
nl = [(1:branch_num)'; (1:branch_num)'];     
Cft = sparse(nl, [f;t], [ones(branch_num, 1); -ones(branch_num, 1)],branch_num , bus_num);
Bf = sparse(nl, [f; t], [b; -b]); %branch_num x bus_num
Bbus = Cft' * Bf; %  bus_num x bus_num
power_gen=mpc.gen(:,1);
Ag=sparse(power_gen,1:gen_num,ones(gen_num,1),bus_num,gen_num);

slack_bus = 1;  
noref   = (2:bus_num)';      
noslack = find((1:bus_num)' ~= slack_bus);

PTDF=makePTDF(mpc.baseMVA,mpc.bus,mpc.branch,4);



%% the tightening of large constants technique
Mbrup=BR-branch_max*ones(1,num);
Mbrdn=-BR-branch_max*ones(1,num);
Mgup=PP-P_max*ones(1,num);
Mgdn=-PP+P_min*ones(1,num);
%% the screening of redundant constraints technique
for br=1:branch_num
for   n=1:num
  if Mbrup(br,n)<0
      Mbrup(br,n)=0;
  end
  if Mbrdn(br,n)<0
      Mbrdn(br,n)=0;
  end
end
end
for i=1:gen_num
for   n=1:num
  if Mgup(i,n)<0
      Mgup(i,n)=0;
  end
  if Mgdn(i,n)<0
      Mgdn(i,n)=0;
  end
end
end


brup=[8    36    38    93   119   126   127   128   129   141   161 163 177]';
brdn=[7  9  79    94    96    97   104   137 177]';
gup=[1 2 6 7  9 11 16 19 21    26  29 30 36 37 38 44 48  49 52 54]';

gdn=[1:54]';

brup=[brup;zeros(branch_num-size(brup,1),1)];
brdn=[brdn;zeros(branch_num-size(brdn,1),1)]; 
gup=[gup;zeros(gen_num-size(gup,1),1)];  
gdn=[gdn;zeros(gen_num-size(gdn,1),1)];  



% ee=sort([8;104]);

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
% Ag_w=sparse(wind_gen,1:wind_num,ones(wind_num,1),bus_num,wind_num);
slack_bus = 1;   %平衡节点设置
noref   = (2:bus_num)';      %将1作为电压参考节点
noslack = find((1:bus_num)' ~= slack_bus);
PTDF= zeros(branch_num, bus_num);
PTDF(:, noslack) = full(Bf(:, noref) / Bbus(noslack, noref));

% Beta=mpc.gen(:,9)/sum(mpc.gen(:,9));
% Beta=[0.1;0.2;0.2;0.4;0.1];
% for i=1:bus_num
%    for j=1:branch_num
%       K(j,i)=-PTDF(j,:)*Ag*Beta+PTDF(j,i); 
%    end
%     
% end

PTDF=makePTDF(mpc.baseMVA,mpc.bus,mpc.branch,4);





%% optimization variable definition
MM=500;
Beta=sdpvar(gen_num,1,'full');
p=sdpvar(gen_num,1,'full');
z_br_up=binvar(num,branch_num,'full');
z_br_dn=binvar(num,branch_num,'full');
z_gen_up=binvar(num,gen_num,'full');
z_gen_dn=binvar(num,gen_num,'full');
for i=1:bus_num
   KA(:,i)=-PTDF*Ag*Beta+PTDF(:,i); 
end
KA=-KA;
%% objective function

fmin=[mpc.gencost(:,5)*100+mpc.gencost(:,6)]'*p;
%% constraints

lv=[sum(PD)-sum(p)==0];  % power balance constraint
branch_cons=[]; % transmission line constraints
gen_cons=[]; % generator output constraints

%% up tranismission line constraints
kk=1;
for  br=1:branch_num
  if br==brup(kk) 
for n=1:num
%     if Mbrup(br,n)>0
   branch_cons=[branch_cons,PTDF(br,:)*(Ag*p-PD)+KA(br,:)*U_bus(:,n)<=branch_max(br)+Mbrup(br,n)*z_br_up(n,br)]; %%
%     end
end
  kk=kk+1;
  else
  branch_cons=[branch_cons,PTDF(br,:)*(Ag*p-PD)<=branch_max(br)] ;
  end
end

kk=1;
for br=1:branch_num
    if br==brup(kk)
   kk=kk+1;
branch_cons=[branch_cons,sum(z_br_up(:,br))<=num*eeee];
    end
end
%% down tranismission line constraints
kk=1;
for  br=1:branch_num
  if br==brdn(kk)
   kk=kk+1;
for n=1:num

 branch_cons=[branch_cons,PTDF(br,:)*(Ag*p-PD)+KA(br,:)*U_bus(:,n)>=-branch_max(br)-Mbrdn(br,n)*z_br_dn(n,br)];

end
  else
  branch_cons=[branch_cons,PTDF(br,:)*(Ag*p-PD)>=-branch_max(br)]; 
  end
end
kk=1;
for br=1:branch_num
     if br==brdn(kk)
   kk=kk+1;
branch_cons=[branch_cons,sum(z_br_dn(:,br))<=num*eeee];
    end 
end
%% max generator output constraints
kk=1;
for  i=1:gen_num
  if i==gup(kk)
   kk=kk+1;
for n=1:num

   gen_cons=[gen_cons,p(i)+Beta(i)*ones(1,bus_num)*U_bus(:,n)<=P_max(i)+Mgup(i,n)*z_gen_up(n,i)];%

end
  else
   gen_cons=[gen_cons,p(i)<=P_max(i)];
  end
end
kk=1;
for i=1:gen_num
    if i==gup(kk)
   kk=kk+1;
gen_cons=[gen_cons,sum(z_gen_up(:,i))<=num*eeee];
    end
end

%% min generator output constraints
kk=1;
for  i=1:gen_num
  if i==gdn(kk)
   kk=kk+1;
for n=1:num

   gen_cons=[gen_cons,p(i)+Beta(i)*ones(1,bus_num)*U_bus(:,n)>=P_min(i)-Mgdn(i,n)*z_gen_dn(n,i)];%

end
  else
    gen_cons=[gen_cons,p(i)>=P_min(i)];%
  end
end
kk=1;
for i=1:gen_num
    if i==gdn(kk)
   kk=kk+1;
gen_cons=[gen_cons,sum(z_gen_dn(:,i))<=num*eeee];
    end
end


b_cons=[];
b_cons=[b_cons,Beta>=0];
b_cons=[b_cons,0.99999<=sum(Beta)<=1.00001];
cons=[lv,branch_cons,gen_cons,b_cons];


ops=sdpsettings('solver', 'gurobi','savesolveroutput',1,'gurobi.timelimit',500,'gurobi.MIPGap',0);
tic
result=optimize(cons,fmin,ops);
solution_time=toc;
Afmin=value(fmin) ;

app=value(p);
p=value(p);
Beta=value(Beta);
%% violation probability caculation

DD=bus_p-bus_w;
AAp=Beta.*sum(DD-PD,1);
p=p+AAp;
k=0;

br=PTDF*(Ag*p-DD);
br=br(:,1:num);

for i=1:size(br,1)
  Ab_up(i)=size(find(br(i,1:num)>branch_max(i)+0.001),2)/num;
  Ab_dn(i)=size(find(br(i,1:num)<-branch_max(i)-0.1),2)/num;
end


for i=1:size(p,1)
Ag_up(i)=size(find(p(i,1:num)>P_max(i)),2)/num;
Ag_dn(i)=size(find(p(i,1:num)<0-0.001),2)/num;

end

Abr=max([Ab_up,Ab_dn]);
Agg=max([Ag_up,Ag_dn]);
z_br_dn=value(z_br_dn);
z_br_up=value(z_br_up);

