%% Run only homoclinic (assumes that run 'hill_longpo_beta=<thisbeta>' exists
%Here we run the hilltop bifurcation problem whilst making the lambda
%parameter be the difference away from the curve of Hopf bifurcations

%original equations are
%x_dot=x^2-lambda+2*alpha*y-gamma
%y_dot=y^2-lambda+2*beta*x+gamma

%curve of Hopf bifurcations/neutral saddles is
%lambda_tr0=gamma*(gamma+alpha^2-beta^2)/(alpha+beta)^2

%we set mu=lambda-lambda_tr0


%------------------------
%% Define the system
%------------------------
hill_top_def;
runidbet=sprintf('hill_longpo_beta=%3.2f',thisbeta);
eqnames={'tr1','tr2','fn1','xeq1','yeq1','t1','det1','fn2','xeq2','yeq2','t2','det2'};
runidhom=sprintf('hill_hom_beta=%3.2f',thisbeta);

%% The family of high-period periodic orbits of constant period
% approximates a family of homoclinic connections to a saddle equilibrium.
% additional variables are stored and recorded along the branch:
% tr1,tr2: transition times, when f_1 is minimal/maximal
% fn1,fn2: minimal norm of f between transition times: the times t1,t2 at which
% this occurs are equilibria
% xeq1,yeq1: approx equilibrium 1 (at time t1)
% det1: determinant of jacobian in (xeq1,yeq1)
% xeq2,yeq2,t2,det2: same for equilibrium 2
prob = coco_prob();
prob = coco_set(prob, 'po', 'bifus', 'off');
prob = ode_po2po(prob, '', runidbet, 2);
prob = coco_xchg_pars(prob, 'gamma', 'po.period');
prob = po_add_func(prob, '', 'saddle_q', ...
  @(data,xbp,T0,T,p)saddle_q(data,xbp,T0,T,p,funcs,1,1e-2), eqnames, 'regular');
prob = coco_set(prob, 'cont', 'NAdapt', 1,'NPR',1,'norm', inf,'h_max',100,'PtMX', [5500,50]);

fprintf('\n Following the curve of homoclinic orbits')

coco(prob, runidhom, [], 1, [{'mu' 'gamma' 'po.period'},eqnames], [-10 10]);
%%
figure(1)
hold on
thm = struct('ustab', '', 'xlab', '\mu', 'ylab', '\gamma');
coco_plot_bd(thm,runidhom,'mu', 'gamma','MAX(x)')
coco_plot_bd(thm,runidhom,'mu', 'gamma','MIN(x)')

