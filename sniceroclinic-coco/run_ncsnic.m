%% Continue non-central SNIC in 3 parameters (assumes that run 'hill_hom_beta=<thisbeta>' exists
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
format compact
addpath([pwd(),'/coll_mesh']);
hill_top_def;
runidhom=sprintf('hill_hom_beta=%3.2f',thisbeta);
hombd=coco_bd_table(runidhom);
sn1=coco_bd_table('hill_sn_run1');
labs=coco_bd_labs(runidhom);
idx=coco_bd_lab2idx(runidhom,labs);
%% pick a point along the curve of long-period approx homoclinics that coincides on SN
% near mu=0.05, gamma=3.6
figure(1);
plot(hombd.mu,hombd.gamma,'.-',sn1.mu,sn1.gamma,'.-');
%%
snic_id=find(diff(sign(hombd.mu+0.05)));
snic_id=snic_id(hombd.gamma(snic_id)>3.5);
snic_pt=hombd(snic_id,:);
snic_lab=snic_pt.LAB;
[psol,u0,cseg]=coll_from_sol('po.orb.coll',runidhom,snic_lab);
[seg,ueq,dev,S,V0,W0]=coll_split(psol,funcs.dfdx);
newnames={'xeq','yeq','s1','s2'};%,'dist'};
varnames=[params,{'x0','y0','x1','y1'},newnames];
ic=[varnames;num2cell(1:length(varnames))];
iv=struct(ic{:});
snic_data=struct('iv',iv,'ip',ip,'funcs',funcs);
y0_bdsnic([iv.mu,iv.gamma,iv.alpha,iv.beta])=psol.parameter(:);
y0_bdsnic=y0_bdsnic(:);
[u0,u1]=deal(seg.profile(:,1),seg.profile(:,end));
y0_bdsnic([  iv.xeq,  iv.yeq])=ueq;
y0_bdsnic([iv.x0,iv.y0, iv.x1,iv.y1,iv.s1,iv.s2])=...
          [u0(1), u0(2),u1(1),u1(2),dev(1),dev(2)];
[~,res]=snic_res('',snic_data,y0_bdsnic);
disp('initial residual at boundary')
disp(res);
tini0=[reshape(psol.mesh(1:end-1,:),[],1);psol.mesh(end)];
tini=tini0*psol.period;
uini=coll_eva(psol,tini0);
isolargs={tini(:),uini',psol.parameter(:)};
shvals=struct('xeq',ueq(1),'yeq',ueq(2),'s1',dev(1),'s2',dev(2));
%%
prob=init_snic(coco_prob(),'init',snic_data,'isolargs',isolargs,'table',shvals);
prob = coco_set(prob, 'cont', 'NAdapt', 1,'NPR',1,'norm', inf,'PtMX', [20,20]);
fprintf('\n correcting snic in 1 parameter\n')
coco(prob, 'snic_ini', [], 0, {'mu','s1','s2','xeq','yeq','x1','y1','T'}, [-10 10]);
%% reload run, now increasing T to reduce abs(s2) freeing T and dropping phase
eplabs=coco_bd_labs('snic_ini','EP');
prob=init_snic(coco_prob(),'reload',snic_data,'run','snic_ini','lab',eplabs(2));
prob = coco_set(prob, 'cont', 'NAdapt', 1,'NPR',100,'norm', inf,'h_max',100,'PtMX', [0,1200]);
fprintf('\n Reduce s2 by increasing T\n')
coco(prob, 'snic_s2', [], 1, {'T','s2','s1','mu','xeq','yeq','x1','y1'}, [170 2000]);
%% reload run, change mu, gamma, detect s2=0
eps2labs=coco_bd_labs('snic_s2','EP');
prob=init_snic(coco_prob(),'reload',snic_data,'run','snic_s2','lab',eps2labs(2));
prob = coco_add_event(prob, 'NCSNIC','boundary', 's2','>',0);
prob = coco_set(prob, 'cont', 'NAdapt', 1,'NPR',1,'norm', inf,'PtMX', [10,20]);
fprintf('\n vary two parameters, check s2\n')
coco(prob, 'snic2nc', [], 1, {'mu','gamma','s1','s2','xeq','yeq','x1','y1','dist'}, [-0.05,0]);
%%
animate_timeprofiles(iv,ip,funcs,'snic2nc',sn1,hombd);
%% reload run, continue NCSNIC in mu, gamma, beta, fixing s2=0
nclabs=coco_bd_labs('snic2nc','NCSNIC');
prob=init_snic(coco_prob(),'reload',snic_data,'run','snic2nc','lab',nclabs);
prob = coco_add_event(prob, 'invalid', 'boundary','dist','>',1e0);
prob = coco_add_event(prob, 'NCSNIC', 'boundary','mu','>',0);
prob = coco_set(prob, 'cont', 'NAdapt', 1,'NPR',1,'norm', inf,'PtMX', 2000,'h0',1e-3);
fprintf('\n vary three parameters\n')
coco(prob, 'ncsnic', [], 1, {'mu','gamma','beta','s1','xeq','yeq','x1','y1','dist'}, {[-1.2,0],[-9,5],[-5,0]});
