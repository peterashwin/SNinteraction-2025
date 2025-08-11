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

syms x y mu alpha beta gamma
f=[x^2-(mu+gamma*(gamma+alpha^2-beta^2)/((alpha+beta)^2))+2*alpha*y-gamma;
    y^2-(mu+gamma*(gamma+alpha^2-beta^2)/((alpha+beta)^2))+2*beta*x+gamma];

F = sco_sym2funcs(f, {[x; y], [mu; gamma; alpha; beta]}, ...
  {'x', 'p'}, 'maxorder', 3, 'filename', 'sys_hill');


%assign the vector field and derivative to funcs
funcs = struct('f', F(''), 'dfdx', F('x'), 'dfdp', F('p'), ...
  'dfdx_dx', F('x*v'), 'dfdxdx_dx', F({'x','x*v'}), ...
  'dfdpdx_dx', F({'p','x*v'}));


%------------------------
%% Find equilibria and continue Hopf bifurcation curve
%------------------------

% continue equilibria to find Hopf and saddle-node bifurcations
fprintf('\n Following a branch of equilibria')

params = {'mu' 'gamma' 'alpha' 'beta'};
%initial parameter values
thisbeta=-1.3;
p0 = [1; 1; 3.1; thisbeta];
%initial guess at equilibrium
x0=[-1 ; 1]; 

prob = coco_prob;
prob = ode_isol2ep(prob, '', funcs, x0, params, p0);
coco(prob, 'hill_ep_run', [], 1, 'mu', [-1 1]);

% continue the Hopf bifurcations in mu-gamma space (should be at mu=0!)
fprintf('\n Following the curve of Hopf bifurcations')

HB   = coco_bd_labs('hill_ep_run', 'HB');
prob = coco_prob;
prob = ode_HB2HB(prob, '', 'hill_ep_run', HB);
data = struct('dfdxhan', F('x'), 'Dfdxdxhan', F({'x*v','x*v'}), ...
  'Dfdxdxdxhan', F({'x*v','x*v','x*v'}), 'nanflag', 1);

prob = coco_add_event(prob, 'UZ', 'gamma',1.0);

coco(prob, 'hill_hb_run', [], 1, {'mu', 'gamma', 'L1'}, [-1 1]);

figure(1)
%clf
hold on
thm = struct('special', {{'BTP'}});
thm.xlab = '\mu';
thm.ylab = '\gamma';
thm.zlab = 'x_1';
coco_plot_bd(thm, 'hill_hb_run', 'mu', 'gamma', 'MAX(x)')
grid on

%---------------------------------------------------------
%% continue saddle-node bifurcations in mu-gamma space
%---------------------------------------------------------
BTP = coco_bd_labs('hill_hb_run', 'BTP');

figure(1)
hold on
thm.special = {};
for i=1:length(BTP)
  prob = coco_prob;
  prob = ode_ep2ep(prob, '', 'hill_hb_run', BTP(i), '-no-var');
  sol  = ep_read_solution('hill_hb_run', BTP(i));
  runid = [sprintf('hill_ep_run_gamma=%3.2f', sol.p(2)),sprintf('_beta=%3.2f',sol.p(4))];
  coco(prob, runid, [], 1, 'mu', [-1 1]);
  SN = coco_bd_labs(runid, 'SN');
  prob = coco_prob;
  prob = ode_SN2SN(prob, '', runid, SN(1));
  runid = sprintf('hill_sn_run%d', i);
  coco(prob, runid, [], 1, {'mu','gamma'}, [-1 1]);
  coco_plot_bd(thm, runid, 'mu', 'gamma', 'MAX(x)')
  drawnow
end
hold off

%---------------------------------------------
%% Continue families of periodic orbits
%---------------------------------------------
UZlabs = coco_bd_labs('hill_hb_run', 'UZ');

figure(1)
hold on
for lab=UZlabs
  prob = coco_prob;
  prob = ode_ep2ep(prob, '', 'hill_hb_run', lab, '-no-var');
  sol = ep_read_solution('hill_hb_run', lab);
  runid = [sprintf('hill_ep_run_gamma=%3.2f', sol.p(2)),sprintf('_beta=%3.2f',sol.p(4))];
  coco(prob, runid, [], 1, 'mu', [-1 1]);
  thm.special = {'HB', 'SN'};
  coco_plot_bd(thm, runid, 'mu', 'gamma', 'x')
  drawnow
  
  chart = coco_read_solution('hill_hb_run', lab, 'chart');
  %if chart.p(6)<0 % supercritical Hopf bifurcations
    prob = coco_prob;
    prob = coco_set(prob, 'po', 'bifus', 'on');
    prob = ode_HB2po(prob, '', 'hill_hb_run', lab);
    sol = ep_read_solution('hill_hb_run', lab);
    runid = [sprintf('hill_po_run_gamma=%3.2f', sol.p(2)),sprintf('_beta=%3.2f',sol.p(4))];
    prob = coco_set(prob, 'cont', 'NAdapt', 1, 'PtMX', [500, 0]);
    coco(prob, runid, [], 1, {'mu' 'po.period'}, ...
      {[-1 1], [0 50]});
    thm.special = {'SN'};
    coco_plot_bd(thm, runid, 'mu', 'gamma', 'MAX(x)')
    coco_plot_bd(thm, runid, 'mu', 'gamma', 'MIN(x)')
    drawnow
  %end
end
hold off

%---------------------------------------------
%% Continue the curve of saddle node of periodic orbits
%---------------------------------------------

%these are the parameters for the po run with the SNPO
gsnpo=1;
betasnpo=thisbeta;
runid=[sprintf('hill_po_run_gamma=%3.2f', gsnpo),sprintf('_beta=%3.2f',thisbeta)];

SN = coco_bd_labs(runid, 'SN');
prob = coco_prob;

prob = ode_SN2SN(prob, '', runid, min(SN));
prob = coco_set(prob, 'cont', 'norm', inf, 'h_max', 200, 'MaxRes', 1, ...
  'NAdapt', 1, 'NPR', 10, 'PtMX', 1000);

coco(prob,'hill_snpo_run', [], 1, {'mu', 'gamma', 'po.period'}, ...
  {[-2 2]});

figure(1)
hold on
coco_plot_bd('hill_snpo_run', 'mu', 'gamma', 'MAX(x)')
hold off


%------------------------------------------------------------------------
%% Continue the curve of homoclinics (approximated by high-period orbits)
%------------------------------------------------------------------------

% Locate a high-period periodic orbit that approximates a homoclinic orbit

%these are the parameters for the po run with the approx hom
gsnpo=1;
betasnpo=thisbeta;
runid=[sprintf('hill_po_run_gamma=%3.2f', gsnpo),sprintf('_beta=%3.2f',thisbeta)];
EP=coco_bd_labs(runid,'EP');

[sol, data] = coll_read_solution('po.orb', runid, max(EP)); %Should be the last PO in the branch. 


f = feval(funcs.f,sol.xbp', repmat(sol.p, [1 size(sol.xbp, 1)])); % Evaluate vector field at basepoints
[~, idx] = min(sqrt(sum(f.*f, 1))); % Find basepoint closest to equilibrium

scale = 2;
%scale=50;
T  = sol.T;
t0 = [sol.tbp(1:idx,1) ; T*(scale-1)+sol.tbp(idx+1:end,1)]; % Crank up period by factor scale
x0 = sol.xbp;
p0 = sol.p;

% Initialize continuation problem structure with the same number of
% intervals as in previous run.

prob = coco_prob();
prob = coco_set(prob, 'coll', 'NTST', data.coll.NTST);
prob = coco_set(prob, 'po', 'bifus', 'off');
prob = ode_isol2po(prob, '', funcs, t0, x0, {'mu' 'gamma' 'alpha' 'beta'}, p0);

prob = coco_set(prob, 'cont', 'NAdapt', 10);
prob = coco_xchg_pars(prob, 'gamma', 'po.period');

runid=sprintf('hill_longpo_beta=%3.2f',thisbeta);

coco(prob, runid, [], 0, {'mu' 'po.orb.coll.err_TF' 'po.period'});

% The family of high-period periodic orbits of constant period
% approximates a family of homoclinic connections to a saddle equilibrium.

prob = coco_prob();
prob = coco_set(prob, 'po', 'bifus', 'off');
prob = ode_po2po(prob, '', runid, 2);
prob = coco_xchg_pars(prob, 'gamma', 'po.period');
prob = coco_set(prob, 'cont', 'NAdapt', 1,'NPR',50,'norm', inf,'h_max',100,'PtMX', 5000);

fprintf('\n Following the curve of homoclinic orbits')

runid=sprintf('hill_hom_beta=%3.2f',thisbeta);
coco(prob, runid, [], 1, {'mu' 'gamma' 'po.period'}, [-10 10]);

figure(1)
hold on
thm = struct('ustab', '', 'xlab', '\mu', 'ylab', '\gamma');
coco_plot_bd(thm,runid,'mu', 'gamma','MAX(x)')
coco_plot_bd(thm,runid,'mu', 'gamma','MIN(x)')




