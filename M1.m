clear
load('case118_uncertaintydata_100.mat');

num=100; % scenario number
%% power system definition

bus_p=bus_p(:,1:num);
bus_w=bus_w(:,1:num);

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

U_bus=bus_p-bus_w-PD; % uncertainty of buses

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

slack_bus = 1;   %
noref   = (2:bus_num)';      
noslack = find((1:bus_num)' ~= slack_bus);


PTDF=makePTDF(mpc.baseMVA,mpc.bus,mpc.branch,4);




%% optimization variable definition
MM=1000;
Beta=sdpvar(gen_num,1,'full');
p=sdpvar(gen_num,1,'full'); %generator outputs
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

lv=[sum(PD)-sum(p)==0]; % power balance constraint
branch_cons=[]; % transmission line constraints
gen_cons=[];% generator output constraints
for n=1:num

for  br=1:branch_num

   branch_cons=[branch_cons,PTDF(br,:)*(Ag*p-PD)+KA(br,:)*U_bus(:,n)>=-branch_max(br)-MM*z_br_dn(n,br)];
   branch_cons=[branch_cons,PTDF(br,:)*(Ag*p-PD)+KA(br,:)*U_bus(:,n)<=branch_max(br)+MM*z_br_up(n,br)]; %%

end
end

for br=1:branch_num
branch_cons=[branch_cons,sum(z_br_up(:,br))<=num*0.05];
branch_cons=[branch_cons,sum(z_br_dn(:,br))<=num*0.05];

end
for i=1:gen_num
  for  n=1:num

gen_cons=[gen_cons,p(i)+Beta(i)*ones(1,bus_num)*U_bus(:,n)<=P_max(i)+MM*z_gen_up(n,i)];
gen_cons=[gen_cons,p(i)+Beta(i)*ones(1,bus_num)*U_bus(:,n)>=P_min(i)-MM*z_gen_dn(n,i)];%
  end
 
end
for i=1:gen_num
gen_cons=[gen_cons,sum(z_gen_up(:,i))<=num*0.05];
gen_cons=[gen_cons,sum(z_gen_dn(:,i))<=num*0.05]; 
end


b_cons=[];
b_cons=[b_cons,Beta>=0];
b_cons=[b_cons,sum(Beta)==1];
cons=[lv,branch_cons,gen_cons,b_cons];


ops=sdpsettings('solver', 'gurobi','savesolveroutput',1,'gurobi.timelimit',500,'gurobi.MIPGap',0);
tic
result=optimize(cons,fmin,ops);
solution_time=toc
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
  Ab_up(i)=size(find(br(i,1:num)>branch_max(i)),2)/num;
  Ab_dn(i)=size(find(br(i,1:num)<-branch_max(i)),2)/num;
end


for i=1:size(p,1)
Ag_up(i)=size(find(p(i,1:num)>P_max(i)),2)/num;
Ag_dn(i)=size(find(p(i,1:num)<0-0.001),2)/num;

end

Abr=max([Ab_up,Ab_dn]);
Agg=max([Ag_up,Ag_dn]);
z_gen_dn=value(z_gen_dn);
z_gen_up=value(z_gen_up);