%% Continue SNICeroclinic in 3 parameters (assumes that run 'hill_hom_beta=<thisbeta>' exists
%original equations are
%x_dot=x^2-lambda+2*alpha*y-gamma
%y_dot=y^2-lambda+2*beta*x+gamma

%curve of Hopf bifurcations/neutral saddles is
%lambda_tr0=gamma*(gamma+alpha^2-beta^2)/(alpha+beta)^2

%we set mu=lambda-lambda_tr0


%------------------------
%% Define the system
%------------------------
clear
hill_top_def;
runidhom=sprintf('hill_hom_beta=%3.2f',thisbeta);
hombd=coco_bd_table(runidhom);
labs=coco_bd_labs(runidhom);
idx=coco_bd_lab2idx(runidhom,labs);
%%
[~,idxmin]=min(abs(abs((hombd{idx,'tr1'}-hombd{idx,'tr2'})')-0.5));
%[fnmin,idxmin]=min(max(hombd{idx,{'fn1','fn2'}},[],2));
imin=idx(idxmin);
pohomsnic=po_read_solution(runidhom,labs(idxmin));
shvals=hombd(imin,:);
[~,i1]=min(abs(shvals.t1-pohomsnic.tbp/pohomsnic.T));
[~,i2]=min(abs(shvals.t2-pohomsnic.tbp/pohomsnic.T));
assert(abs(shvals.det1),abs(shvals.det2)); % confirm that eq1 is SN, eq2 is saddle
newnames={'xeq1','yeq1','xeq2','yeq2','s1','s2'};
varnames=[params,{'x0','y0','x1','y1'},newnames];
ic=[varnames;num2cell(1:length(varnames))];
iv=struct(ic{:});
xsn=[shvals.xeq1,shvals.yeq1];
xsa=[shvals.xeq2,shvals.yeq2];
irg=i1+1:i2-1;
uini=pohomsnic.xbp(irg,:);
tini=pohomsnic.tbp(irg);
tini=tini-tini(1);
pini=pohomsnic.p;
T=tini(end);
homsnic_data=struct('iv',iv,'ip',ip,'funcs',funcs);
isolargs={tini,uini,params,pini};
%%
figure(2);clf;
plot(tini,uini,'o-');
prob=init_homsnic(coco_prob(),'init',homsnic_data,'isolargs',isolargs,'table',shvals);
prob=add_hetphasecond(prob,'seg.coll');
prob = coco_set(prob, 'cont', 'NAdapt', 1,'NPR',1,'norm', inf,'h_max',100,'PtMX', [20,20]);
%%
fprintf('\n correcting homsnic in 2 parameters\n')
coco(prob, 'homsnic_ini', [], 0, [{'mu', 'gamma'},newnames], [-10 10]);
%% reload run, now increasing T to reduce s1,s2 freeing T
eplabs=coco_bd_labs('homsnic_ini','EP');
prob=init_homsnic(coco_prob(),'reload',homsnic_data,'run','homsnic_ini','lab',eplabs(2));
prob=add_hetphasecond(prob,'seg.coll');
prob = coco_set(prob, 'cont', 'NAdapt', 1,'NPR',5,'norm', inf,'h_max',100,'PtMX', [0,50]);
fprintf('\n Reduce s1 and s2 by increasing T\n')
coco(prob, 'homsnic_s2', [], 1, [{'T','mu', 'gamma'},newnames], [60 1000]);
%% reload run, continue homsnic in 3 parameters fixing T, phase, leaving s1,s2 free
eplabs=coco_bd_labs('homsnic_s2','EP');
prob=init_homsnic(coco_prob(),'reload',homsnic_data,'run','homsnic_s2','lab',eplabs(2));
prob=add_hetphasecond(prob,'seg.coll');
prob = coco_set(prob, 'cont', 'NAdapt', 1,'NPR',1,'norm', inf,'h_max',10,'PtMX', [-30,50]);
fprintf('\n Continue homsnic in 3 parameters`n')
coco(prob, 'homsnic', [], 1, [{'mu', 'gamma','beta'},newnames], [-20,0]);
%% check return trajectory
homsnic=coco_bd_table('homsnic');
npt=size(homsnic,1);
irg=1:npt;
T=100;
tret=linspace(0,T,100*T+1);
clear t_return u_return
for i=1:length(irg)
    lab=homsnic.LAB(irg(i));
    fprintf('extracting i=%d of %d,lab=%d\n',i,npt,lab);
    ch=coco_read_solution('homsnic','homsnic',lab,'chart');
    p=ch.x(ismember(fieldnames(iv),params));
    [usn,usa]=deal([ch.x(iv.xeq1);ch.x(iv.yeq1)],[ch.x(iv.xeq2);ch.x(iv.yeq2)]);
    [vsa,dsa]=eig(funcs.dfdx(usa,p));
    dsa=diag(dsa);
    vsa=vsa(:,dsa>0);
    vsa=-vsa*sign(vsa(1));
    sol=ode45(@(t,x)funcs.f(x,p),[0,T],usa+vsa*1e-1,odeset('RelTol',1e-8,'AbsTol',1e-8));
    uret=deval(sol,tret);
    t_return{irg(i)}=tret;
    u_return{irg(i)}=uret;
end
homsnic.('t_return')=t_return(:);
homsnic.('u_return')=u_return(:);
save('homsnic.mat','homsnic');