%% Plot bifurcation diagram (Figure 2)
clear 
format compact
startup_coco(fullfile('coco_r3316'))
exportweb=true;
homsnic=coco_bd_table('homsnic');
ncsnic=coco_bd_table('ncsnic');
btp=coco_bd_table('hill_btp_run');
gh=coco_bd_table('hill_gh_run');
%% cusp formulas from paper
lamtr0=@(gamma,beta,alpha)gamma.*(gamma+alpha.^2-beta.^2)./((alpha+beta).^2);
lamcuf=@(alpha,beta)3/2*alpha.^(2/3).*(beta.^2).^(1/3).*(alpha.^(2/3)+(beta.^2).^(1/3));
gamcuf=@(alpha,beta)3/2*alpha.^(2/3).*(beta.^2).^(1/3).*(alpha.^(2/3)-(beta.^2).^(1/3));
cu.beta=linspace(0,-3.1,100);
cu.alpha=gh.alpha(1)+0*cu.beta;
cu.gamma=gamcuf(cu.alpha,cu.beta);
cu.lambda=lamcuf(cu.alpha,cu.beta);
cu.mu=cu.lambda-lamtr0(cu.gamma,cu.beta,cu.alpha);
%% plot
clr=lines();
c=struct('TB','k','gh',clr(3,:),'sniceroclinic',clr(1,:),'cusp',clr(2,:),...
    'ncsnic',clr(5,:));
[lw,lw3,lw4]=deal({'LineWidth',2},{'LineWidth',3},{'LineWidth',4});
txt={'FontSize',16,'FontName','Courier'};
ltx={'Interpreter','LaTeX'};
fig=figure(1);clf;ax=gca;hold(ax,'on');
[lims.x,lims.y,lims.z]=deal([-0.12,0.04],[-3.5,0],[-4,4]);
plot3(ax,homsnic.mu,homsnic.beta,homsnic.gamma,'color',c.sniceroclinic,...
    'DisplayName','SNICeroclinic',lw{:});
plot3(ax,ncsnic.mu,ncsnic.beta,ncsnic.gamma,'color',c.ncsnic,...
    'DisplayName','non-central SNIC',lw{:});
plot3(ax,btp.mu,btp.beta,btp.gamma,'color',c.TB,...
    'DisplayName','Takens-Bogdanov bifurcation',lw{:});
plot3(ax,gh.mu,gh.beta,gh.gamma,'color',c.gh,...
    'DisplayName','Bautin bifurcation',lw{:});
plot3(ax,cu.mu,cu.beta,cu.gamma,'color',c.cusp,...
    'DisplayName','Cusp bifurcation',lw{:});
grid(ax,'on');
set(ax,txt{:});
xlim(ax,lims.x);
ylim(ax,lims.y);
zlim(ax,lims.z);
ax.View=[20,30];
xlabel(ax,'$\mu$',txt{:},ltx{:});
ylabel(ax,'$\beta$',txt{:},ltx{:});
zlabel(ax,'$\gamma$',txt{:},ltx{:});
lgd = legend(ax,'FontName','Times','Interpreter','latex','Location','EastOutside');
lgd.EdgeColor='none';
fig.Position(3:4)=[900,750];
if exportweb && ~verLessThan('matlab','26') %#ok<*VERLESSMATLAB>
    exportgraphics(figure(1),'../Figure5.html');
end
