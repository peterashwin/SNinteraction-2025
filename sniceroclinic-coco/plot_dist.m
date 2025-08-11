%% plot distance of orbit segment ens to equilibrium
clear
homsnic=coco_bd_table('homsnic');
ncsnic=coco_bd_table('ncsnic');
figure(1);clf;ax1=gca;
lw={'linewidth',2};
ltx={'Interpreter','latex','FontSize',20};
txt={'FontSize',20,'FontName','Courier','FontWeight','bold'};
semilogy(homsnic.beta,sum(homsnic{:,{'s1','s2'}},2),'+-',...
    'DisplayName','$|s_-|+|s_+|$ for SNICeroclinic',lw{:});
hold(ax1,'on');
semilogy(ncsnic.beta,ncsnic{:,{'dist'}},'+-',...
    'DisplayName','$\|u_+-u_\mathrm{sn}\|+|s_-|$ for non-central SNIC',lw{:});
legend(ax1,'Location','best',ltx{:});
xlabel(ax1,'$\beta$',ltx{:});
set(ax1,lw{:},txt{:});
grid(ax1,'on');
xlim(ax1,[-3.2,0.1]);
ylim(ax1,[1e-3,1e-1]);
%%
exportgraphics(figure(1),'fig-snics_dist.pdf');