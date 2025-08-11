%% plot homsnics along branch
clear
exportvideo=false;
exportplot=false;
hill_top_def;
s=load('homsnic.mat');
homsnic=s.homsnic;
runidhom=sprintf('hill_hom_beta=%3.2f',thisbeta);
hb=coco_bd_table('hill_hb_run');
snpo=coco_bd_table('hill_snpo_run');
sn1=coco_bd_table('hill_sn_run1');
sn2=coco_bd_table('hill_sn_run2');
hombd=coco_bd_table(runidhom);
homsnicnames={'xeq1','yeq1','xeq2','yeq2','s1','s2'};
homsnicvarnames=[params,{'x0','y0','x1','y1'},homsnicnames];
ic_hs=[homsnicvarnames;num2cell(1:length(homsnicvarnames))];
iv_hs=struct(ic_hs{:});
ncsnic=coco_bd_table('ncsnic');
ncsnicnames={'xeq','yeq','s1','s2','dist'};
ncsnicvarnames=[params,{'x0','y0','x1','y1'},ncsnicnames];
ic_nc=[ncsnicvarnames;num2cell(1:length(ncsnicvarnames))];
iv_nc=struct(ic_nc{:});
if exportvideo
    vid=VideoWriter('homsnic','Motion JPEG AVI');
    vid.FrameRate=5;
    open(vid);
end
figure(2);clf;
tl=tiledlayout(1,2);nexttile;ax1=gca;
lw={'linewidth',2};
ltx={'Interpreter','latex','FontSize',20};
txt={'FontSize',20,'FontName','Courier','FontWeight','bold'};
hold(ax1,'on');
clr=lines();
[s2bd,retbd]=deal(0.2,0.1);
s2large=homsnic.s2>s2bd;
yrdist=cellfun(@(y,xsn,ysn)norm(y(:,end)-[xsn;ysn]),...
    homsnic.u_return,num2cell(homsnic.xeq1),num2cell(homsnic.yeq1));
retdistlarge=yrdist>retbd;
retonly=retdistlarge&~s2large;
ok=~s2large&~retdistlarge;
plot3(ax1,homsnic.mu(retonly),homsnic.gamma(retonly),homsnic.beta(retonly),...
    'o','color',(clr(2,:)+1)/2,lw{:},'DisplayName',...
    sprintf('$\\|u_\\mathrm{ret}(T)-u_\\mathrm{SN}\\|>%g$\n(no cycle)',retbd));
plot3(ax1,homsnic.mu(s2large),homsnic.gamma(s2large),homsnic.beta(s2large),...
    'o','color',(clr(1,:)+1)/2,lw{:},'DisplayName',sprintf('$s_2>%g$',s2bd));
plot3(ax1,homsnic.mu(ok),homsnic.gamma(ok),homsnic.beta(ok),...
    'o-','color',clr(1,:),lw{:},'DisplayName','homsnic');
legend(ax1,ltx{:});
grid(ax1,'on');
view(ax1,[-75,10]);
xlabel(ax1,'$\mu$',ltx{:});
ylabel(ax1,'$\gamma$',ltx{:});
zlabel(ax1,'$\beta$',ltx{:});
zlim(ax1,[-3.2,0.5]);
set(ax1,txt{:},lw{:})
%%
nexttile;ax2=gca;
npt=size(homsnic,1);
if exportplot
    irg=40;
else
    irg=1:npt;
end
pl1init=false;
ok_only=false;
for i=1:length(irg)
    if ~ok(i) && ok_only
        continue
    end
    lab=homsnic.LAB(irg(i));
    [cs,ds]=coll_read_solution('seg','homsnic',lab,'chart','data');
    t=cs.tbp;
    u=cs.xbp;
    ch=coco_read_solution('homsnic','homsnic',lab,'chart');
    plot(ax2,u(:,1),u(:,2),lw{:});
    hold(ax2,'on');
    [usn,usa]=deal([ch.x(iv_hs.xeq1);ch.x(iv_hs.yeq1)],[ch.x(iv_hs.xeq2);ch.x(iv_hs.yeq2)]);
    [vsa,dsa]=eig(funcs.dfdx(usa,cs.p));
    dsa=diag(dsa);
    vsa=vsa(:,dsa>0);
    vsa=-vsa*sign(vsa(1));
    T=100;
    sol=ode45(@(t,x)funcs.f(x,cs.p),[0,T],usa+vsa*1e-1,odeset('RelTol',1e-8,'AbsTol',1e-8));
    yr=deval(sol,linspace(0,T,100*T+1));
    plot(ax2,yr(1,:),yr(2,:),'k-',lw{:});
    plot(ax2,usn(1),usn(2),'ko','MarkerSize',10,'MarkerFaceColor','y');
    plot(ax2,usa(1),usa(2),'rx','MarkerSize',10,'Linewidth',3);
    hold(ax2,'off');
    %xlim(ax2,[-12,12]);
    %ylim(ax2,[-15,10]);
    xlabel(ax2,'$x$',ltx{:});
    ylabel(ax2,'$y$',ltx{:});
    set(ax2, txt{:},lw{:})
    ax2.YLimitMethod='padded';
    ax2.XLimitMethod='padded';    
    pcoords={cs.p(ip.mu),cs.p(ip.gamma),cs.p(ip.beta)};
    title(ax2,...
        sprintf('i=%3d, $\\mu=$%6.3f, $\\gamma=$%6.3f, $\\beta=$%6.3f, $s_1=$%5.3e, $s_2$=%5.3e',...
        i,pcoords{:},ch.x(iv_hs.s1),ch.x(iv_hs.s2)),ltx{:});
    if ~pl1init
        pl1init=true;
        plseg=plot3(ax1,pcoords{:},'ko','MarkerSize',12,'MarkerFaceColor','r','HandleVisibility','off');
    else
        [plseg.XData,plseg.YData,plseg.ZData]=deal(pcoords{:});
    end
    %pause(0.2);
    drawnow
    if exportvideo
        frame=getframe(figure(2));
        writeVideo(vid,frame);
    end
end
if exportvideo
    close(vid);
end
if exportplot
    exportgraphics(figure(2),'homsnic.pdf');
end
