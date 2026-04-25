%% Plot bifurcation diagram (Figure 2)
clear 
format compact
startup_coco(fullfile('coco_r3316'))
exportweb=true;
clr=lines();
c=struct('hopf',2,'sne',5,'snpo',4,'hom',6);
hb=coco_bd_table('hill_hb_run');
ishopf=hb.('ep.test.BTP')>-1e-3;
isbtp=strcmp(hb.TYPE,'BTP');
isgh=strcmp(hb.TYPE,'GH');
sn1=coco_bd_table('hill_sn_run1');
sn2=coco_bd_table('hill_sn_run2');
snpo=coco_bd_table('hill_snpo_run');
hom=coco_bd_table('hill_hom_beta=-1.30');
iscu=strcmp(sn1.TYPE,'FP');
i_sniceroclinic=find(~isnan(hom.det1)&abs(hom.det1)<3e-5,1,'last');
is_snic1=~isnan(hom.det1)&abs(hom.det1)<2e-2;
i_ncsnic_cand=find(diff(is_snic1));
i_ncsnic=i_ncsnic_cand(2);
[~,i_snpo_end]=min(sum([hom.mu,hom.gamma]-[snpo.mu(end),snpo.gamma(end)],2).^2);
[i_csnpo(1),i_csnpo(2)]=deal(min(i_snpo_end,i_sniceroclinic),max(i_snpo_end,i_sniceroclinic));
rg_csnpo=i_csnpo(1):i_csnpo(2);
[lw,lw3,lw4]=deal({'LineWidth',2},{'LineWidth',3},{'LineWidth',4});
txt={'FontSize',16,'FontName','Courier'};
ltx={'Interpreter','LaTeX'};
fig=figure(1);clf;tiledlayout(2,3,'TileSpacing','tight');
tiles={{1,[2,2]},{3,[1,1]}};
[lims.x,lims.y]=deal({[-0.7,0.4],[-0.15,0.05]},{[-6,5],[3.5,3.7]});
for i=1:length(tiles)
    ax(i)=nexttile(tiles{i}{:});hold(ax(i),'on');
    plot(ax(i),hb.mu(ishopf),hb.gamma(ishopf),'-','Color',clr(c.hopf,:),'DisplayName','Hopf bifurcation',lw{:});
    plot(ax(i),hb.mu(isbtp),hb.gamma(isbtp),'kp','DisplayName','Takens-Bogdanov point',lw{:});
    plot(ax(i),hb.mu(isgh),hb.gamma(isgh),'ks','DisplayName','Bautin point',lw{:},'MarkerFaceColor','k');
    plot(ax(i),sn1.mu,sn1.gamma,'-','Color',clr(c.sne,:),'DisplayName','Saddle-node of equilibria',lw{:});
    plot(ax(i),real(sn2.mu),real(sn2.gamma),'-','Color',clr(c.sne,:),'HandleVisibility','off',lw{:});
    plot(ax(i),snpo.mu,snpo.gamma,'-','Color',clr(c.snpo,:),'DisplayName','Saddle-node of periodic orbits (SNPO)',lw4{:});
    plot(ax(i),hom.mu(rg_csnpo),hom.gamma(rg_csnpo),':','Color',clr(c.snpo,:),'DisplayName','SNPO (conjectured)',lw4{:});
    plot(ax(i),hom.mu,hom.gamma,'-','Color',clr(c.hom,:),'DisplayName','Long-period periodic orbit',lw{:});
    plot(ax(i),hom.mu(i_sniceroclinic),hom.gamma(i_sniceroclinic),'k^','DisplayName','SNICeroclinic',lw{:},'MarkerFaceColor','k');
    plot(ax(i),hom.mu(i_ncsnic),hom.gamma(i_ncsnic),'kd','DisplayName','non-central SNIC',lw{:},'MarkerFaceColor','k');
    plot(ax(i),sn1.mu(iscu),sn1.gamma(iscu),'kv','DisplayName','cusp',lw{:},'MarkerFaceColor','k');
    grid(ax(i),'on');
    set(ax(i),txt{:});
    xlim(ax(i),lims.x{i});
    ylim(ax(i),lims.y{i});
    xlabel(ax(i),'$\mu$',txt{:},ltx{:});
    ylabel(ax(i),'$\gamma$',txt{:},ltx{:});
end
lgd = legend(ax(1),'FontName','Times','Interpreter','latex');
lgd.Layout.Tile = 6;   % place legend in the 4th tile of the grid
lgd.EdgeColor='none';
fig.Position(3:4)=[1400,700];
if exportweb && ~verLessThan('matlab','26') %#ok<*VERLESSMATLAB>
    exportgraphics(figure(1),'../Figure2.html');
end
