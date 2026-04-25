%% Continuation of local codimension-two bifurcations
%------------------------
%% Define the system
%------------------------
clear 
format compact
startup_coco(fullfile('coco_r3316'))
hill_top_def;
%% Continue Takens-Bogdanov point in 3 parameters
TB   = coco_bd_labs('hill_hb_run', 'BTP');
prob = coco_prob;
prob = ode_HB2HB(prob, '', 'hill_hb_run', TB(1));
data = struct('dfdx', F('x'), 'Dfdxdx', F({'x*v','x*v'}), ...
    'Dfdxdxdx', F({'x*v','x*v','x*v'}), 'nanflag', 1);
prob = ep_HB_add_func(prob, '', 'btp', @(d,x,p,v,k)k, struct(),'zero');
coco(prob, 'hill_btp_run', [], 1, {'mu', 'gamma', 'beta'}, {[-1 1],[-4,4],[-3.101,0]});
%% Continue GH in 3 parameters
GH   = coco_bd_labs('hill_hb_run', 'GH');
prob = coco_prob;
prob = ode_HB2HB(prob, '', 'hill_hb_run', GH);
data = struct('dfdx', F('x'), 'Dfdxdx', F({'x*v','x*v'}), ...
    'Dfdxdxdx', F({'x*v','x*v','x*v'}), 'nanflag', 1);
prob = ep_HB_add_func(prob, '', 'GH', @lyapunov,data,'zero');
coco(prob, 'hill_gh_run', [], 1, {'mu', 'gamma', 'beta'}, {[-1 1],[-4,4],[-3.101,0]});
