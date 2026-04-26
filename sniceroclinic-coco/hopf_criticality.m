%% test script to compute Hopf criticality using symbolic toolbox
%% initialize coco
%%
clear
startup_coco(fullfile('coco_r3316'))
%------------------------
%% Define the system
%------------------------
params = {'mu', 'gamma', 'alpha', 'beta'};
ic=[params;num2cell(1:length(params))];
ip=struct(ic{:});
syms x y u v
syms(params)
assumeAlso([alpha+beta,alpha,-beta,gamma]>0)
mu=0;
lambda=gamma*(gamma+alpha^2-beta^2)/((alpha+beta)^2);
f=[x^2-mu+lambda+2*alpha*y-gamma;
        y^2-mu+lambda+2*beta*x+gamma];
%% insert equilibria for Hopf bifurcation
xeq=-gamma/(alpha+beta);
yeq=-xeq;
%% transform coordinates to normal form at first order
Axy=simplify(subs(jacobian(f,[x;y]),[x;y],[xeq;yeq]));
om=sqrt(det(Axy));
ev=sym([1;0]);
V=[om*ev,Axy*ev];
nV2=sum(V.^2,1);
uveq=V\[xeq;yeq];
fuv=simplify(V\subs(f,[x;y],V*[u;v]));
[P,Q]=deal(fuv(1),fuv(2));
[         Puu,         Puv,        Pvv,        Puuu,         Puvv]=deal(...
    diff(P,u,u),diff(P,u,v),diff(P,v,v),diff(P,u,u,u),diff(P,u,v,v));
[          Quu,        Quv,       Qvv,         Quuv,         Qvvv]=deal(...
    diff(Q,u,u),diff(Q,u,v),diff(Q,v,v),diff(Q,u,u,v),diff(Q,v,v,v));
A=simplify(subs(jacobian(fuv,[u;v]),[u;v],uveq));
%% l1 formula from scholarpedia
l1_3=Puuu+Puvv+Quuv+Qvvv;
l1_2=Puv*(Puu+Pvv)-Quv*(Quu+Qvv)-Puu*Quu+Pvv*Qvv;
l1=simplify((l1_3*om+l1_2)/(8*om^2));
%% compare w numerical results
hb=coco_bd_table('hill_hb_run');
hb=hb(hb.('ep.test.BTP')>0,:);
num=@(expr)arrayfun(@(g)double(subs(expr,[mu,alpha,beta,gamma],[0,hb.alpha(1),hb.beta(1),g])),hb.gamma);
l1n=num(l1);   % analytical L1 as a function of gamma
nV2n=num(mean(nV2)); % square of norms of columns in transformation matrix 
fprintf('difference between numerical and analyitcal formula=%g\n',norm(l1n-hb.L1.*nV2n));
%% plot analytical L1 and compare to numerical results
exportweb=true;
clr=lines();
[lw,lw3,lw4]=deal({'LineWidth',2},{'LineWidth',3},{'LineWidth',4});
txt={'FontSize',14,'FontName','Courier'};
ltx={'Interpreter','LaTeX'};
fig=figure(4);clf;tiledlayout(2,1,'TileSpacing','tight');
ax1=nexttile(1);hold(ax1,'on');
plot(ax1,hb.gamma,l1n,'-','Color',clr(1,:),lw{:},'DisplayName','L1 (analytic formula)');
xline(ax1,0,'k','LineWidth',1);
yline(ax1,0,'k','LineWidth',1);
grid(ax1,'on')
ylim(ax1,1*[-1,1]);
L1str=['$$\mbox{Lyapunov coefficient\quad}L_1=',latex(l1),...
    '\quad\mbox{($\alpha=',num2str(hb.alpha(1)),'$, $\beta=',num2str(hb.beta(1)),'$)}$$'];
titlestr=sprintf(['\\textbf{Ashwin et al}, \\textit{Local interaction of two systems with saddle-node '...
    'bifurcations}:\n \\textit{mutualistic and mixed cases},\n %s'],L1str);
title(ax1,titlestr,ltx{:},txt{:});
set(ax1,txt{:});
ax2=nexttile(2);
semilogy(ax2,hb.gamma,abs(l1n-hb.L1.*nV2n),'o','Color',clr(1,:),lw{:});
xlabel(ax2,'$\gamma$',txt{:},ltx{:});
grid(ax2,'on')
title(ax2,'difference to numerical result (after scaling change-of-basis matrix)',ltx{:},txt{:});
ylim(ax2,[1e-18,1]);
set(ax2,txt{:});
fig.Position(3:4)=[900,600];
if exportweb && ~verLessThan('matlab','26') %#ok<*VERLESSMATLAB>
    exportgraphics(figure(4),'../Figure2-Hopf-L1.html');
end


