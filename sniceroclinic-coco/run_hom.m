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

%% animate time profiles
hb=coco_bd_table('hill_hb_run');
snpo=coco_bd_table('hill_snpo_run');
sn1=coco_bd_table('hill_sn_run1');
sn2=coco_bd_table('hill_sn_run2');
hombd=coco_bd_table(runidhom);
labs=coco_bd_labs(runidhom);
idx=coco_bd_lab2idx(runidhom,labs);
clr=lines();
figure(3);clf;tiledlayout(1,3);
txt={'FontSize',18};
ms={'MarkerSize',12};
for i=1:length(labs)
    sol=po_read_solution(runidhom,labs(i));
    bdvals=hombd(idx(i),:);
    nexttile(1);
    ax3=gca;
    plot(ax3,sol.tbp/sol.T,sol.xbp(:,1),'.-','DisplayName','x');
    hold(ax3,'on');
    plot(ax3,sol.tbp/sol.T,sol.xbp(:,2),'.-','DisplayName','y');
    xline(ax3,[bdvals.tr1,bdvals.tr2])
    plot(ax3,bdvals.t1,[bdvals.xeq1,bdvals.yeq1],'ko','MarkerFaceColor',clr(3,:),ms{:});
    plot(ax3,bdvals.t2,[bdvals.xeq2,bdvals.yeq2],'kp','MarkerFaceColor',clr(4,:),ms{:});
    xlabel(ax3,'time');
    ylabel(ax3,'x,y');
    xlim(ax3,[0,1]);
    ylim(ax3,[-8,8]);
    title(ax3,sprintf('mu=%5.2f, gamma=%5.2f',bdvals.mu,bdvals.gamma));
    hold(ax3,'off');
    nexttile(2);ax4=gca;
    plot(ax4,sol.xbp(:,1),sol.xbp(:,2),'.-');
    hold(ax4,'on');
    plot(ax4,bdvals.xeq1,bdvals.yeq1,'ko','MarkerFaceColor',clr(3,:),ms{:});
    plot(ax4,bdvals.xeq2,bdvals.yeq2,'kp','MarkerFaceColor',clr(4,:),ms{:});
    xlabel(ax4,'x');
    ylabel(ax4,'y')
    hold(ax4,'off');
    xlim(ax4,[-8,8]);
    ylim(ax4,[-6,6]);
    nexttile(3);ax5=gca;hold(ax5,'off');
    plot(ax5, sn1.('mu'),sn1.('gamma'),'-','color',clr(1,:),'LineWidth',3);
    hold(ax5,'on');
    %plot(ax5, sn2.('mu'),sn2.('gamma'),'-','color',clr(1,:),'LineWidth',3);
    plot(ax5, hb.('mu'),hb.('gamma'),'r-','LineWidth',3);
    plot(ax5, hombd.('mu'),hombd.('gamma'),'k.-','LineWidth',1);
    plot(ax5,sol.p(ip.mu),sol.p(ip.gamma),'ko','MarkerFaceColor',clr(3,:),ms{:});
    ax5.XLim=[-0.03,0.01];
    ax5.YLim=[3.58,3.61];
    drawnow
end

