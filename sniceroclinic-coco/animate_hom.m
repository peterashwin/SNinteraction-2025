%% after running run_hom, animate time profiles of large-period orbits
clear
hill_top_def;
runidbet=sprintf('hill_longpo_beta=%3.2f',thisbeta);
runidhom=sprintf('hill_hom_beta=%3.2f',thisbeta);
hb=coco_bd_table('hill_hb_run');
snpo=coco_bd_table('hill_snpo_run');
sn1=coco_bd_table('hill_sn_run1');
sn2=coco_bd_table('hill_sn_run2');
hombd=coco_bd_table(runidhom);
labs=coco_bd_labs(runidhom);
idx=coco_bd_lab2idx(runidhom,labs);
clr=lines();
figure(3);clf;tl=tiledlayout(1,3);
for i=1:3
    ax(i)=nexttile(tl,i);
end
[ax1,ax2,ax3]=deal(ax(1),ax(2),ax(3));
txt={'FontSize',18};
ltx={'interpreter','latex'};
ms={'MarkerSize',12};
ax3.XLim=[-0.15,0.05];
ax3.YLim=[3.5,3.75];
for i=1:5:length(labs)
    sol=po_read_solution(runidhom,labs(i));
    bdvals=hombd(idx(i),:);
    plot(ax1,sol.tbp/sol.T,sol.xbp(:,1),'.-','DisplayName','x');
    hold(ax1,'on');
    plot(ax1,sol.tbp/sol.T,sol.xbp(:,2),'.-','DisplayName','y');
    xline(ax1,[bdvals.tr1,bdvals.tr2])
    plot(ax1,bdvals.t1,[bdvals.xeq1,bdvals.yeq1],'ko','MarkerFaceColor',clr(3,:),ms{:});
    plot(ax1,bdvals.t2,[bdvals.xeq2,bdvals.yeq2],'kp','MarkerFaceColor',clr(4,:),ms{:});
    xlabel(ax1,'time');
    ylabel(ax1,'x,y');
    xlim(ax1,[0,1]);
    ylim(ax1,[-8,8]);
    title(ax1,sprintf('mu=%5.2f, gamma=%5.2f, i=%4d',bdvals.mu,bdvals.gamma,i));
    hold(ax1,'off');
    plot(ax2,sol.xbp(:,1),sol.xbp(:,2),'.-');
    hold(ax2,'on');
    plot(ax2,bdvals.xeq1,bdvals.yeq1,'ko','MarkerFaceColor',clr(3,:),ms{:});
    plot(ax2,bdvals.xeq2,bdvals.yeq2,'kp','MarkerFaceColor',clr(4,:),ms{:});
    xlabel(ax2,'x');
    ylabel(ax2,'y')
    hold(ax2,'off');
    xlim(ax2,[-8,8]);
    ylim(ax2,[-6,6]);
    if i==1
        hold(ax3,'on');
        plot(ax3, hb.('mu'),hb.('gamma'),'r-','LineWidth',3);
        plot(ax3, sn1.('mu'),sn1.('gamma'),'-','color',clr(1,:),'LineWidth',3);
        plot(ax3, hombd.('mu'),hombd.('gamma'),'k.-','LineWidth',1);
        po=plot(ax3,sol.p(ip.mu),sol.p(ip.gamma),'ko','MarkerFaceColor',clr(3,:),ms{:});
        xlabel(ax3,'$\mu$',txt{:},ltx{:});
        ylabel(ax3,'$\gamma$',txt{:},ltx{:});        
    else
        [mu,gamma]=deal(sol.p(ip.mu),sol.p(ip.gamma));
        [po.XData,po.YData]=deal(mu,gamma);
        ax3.XLim=[min(ax3.XLim(1),mu-0.05),max(ax3.XLim(2),mu+0.05)];
        ax3.YLim=[min(ax3.YLim(1),gamma-0.2),max(ax3.YLim(2),gamma+0.2)];
    end
    drawnow
end

