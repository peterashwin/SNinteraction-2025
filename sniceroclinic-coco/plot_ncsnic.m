%% plot homsnics along branch
clear
exportvideo=false;
exportplot=false;
hill_top_def;
s=load('homsnic.mat');
homsnic=s.homsnic;
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
    vid=VideoWriter('ncsnic','Motion JPEG AVI');
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
ncok=abs(ncsnic.dist)<1e-1&ncsnic.mu<=0;
plot3(ax1,ncsnic.mu(ncok),ncsnic.gamma(ncok),ncsnic.beta(ncok),...
    '-','color',(clr(4,:)+1)/2,lw{:},'DisplayName',...
    sprintf('NC SNIC'));
plot3(ax1,homsnic.mu(ok),homsnic.gamma(ok),homsnic.beta(ok),...
    '-','color',clr(1,:),lw{:},'DisplayName','homsnic');
legend(ax1,ltx{:},'location','best');
grid(ax1,'on');
view(ax1,[10,10]);
xlabel(ax1,'$\mu$',ltx{:});
ylabel(ax1,'$\gamma$',ltx{:});
zlabel(ax1,'$\beta$',ltx{:});
zlim(ax1,[-3.2,0.5]);
xlim(ax1,[-0.12,0.01]);
set(ax1,txt{:},lw{:})
%%
nexttile;ax2=gca;
%%
labrows=arrayfun(@(i)~isempty(ncsnic{i,'LAB'})&ncok(i),1:length(ncok));
labs=ncsnic{labrows,'LAB'};
npt=length(labs);
if exportplot
    irg=50;
else
    irg=1:npt;
end
pl1init=false;
ok_only=false;
h0=1e-1;
for i=irg
    lab=labs(i);
    [cs,ds]=coll_read_solution('seg','ncsnic',lab,'chart','data');
    t=cs.tbp;
    u=cs.xbp;
    ch=coco_read_solution('snic','ncsnic',lab,'chart');
    plot(ax2,u(:,1),u(:,2),lw{:});
    hold(ax2,'on');
    usn=ch.x([iv_nc.xeq,iv_nc.yeq]);
    [vsn,~,~,J]=snic_eigspace(funcs.dfdx,usn,cs.p,ch.x([iv_nc.x0,iv_nc.y0]));
    ev=eig(J);
    [~,ix]=sort(abs(ev),'descend');
    ev=ev(1);
    plot(ax2,usn(1),usn(2),'ko','MarkerSize',10,'MarkerFaceColor','y');
    plot(ax2,usn(1)+h0*vsn(1)*[-1,1],usn(2)+h0*vsn(2)*[-1,1],'k-',lw{:});
    hold(ax2,'off');
    %xlim(ax2,[-12,12]);
    %ylim(ax2,[-15,10]);
    xlabel(ax2,'$x$',ltx{:});
    ylabel(ax2,'$y$',ltx{:});
    axis(ax2,'equal');
    set(ax2, txt{:},lw{:})
    ax2.YLimitMethod='padded';
    ax2.XLimitMethod='padded';    
    pcoords={cs.p(ip.mu),cs.p(ip.gamma),cs.p(ip.beta)};
    title(ax2,...
        sprintf('i=%3d, $\\mu=$%6.3f, $\\gamma=$%6.3f, $\\beta=$%6.3f, tr ev$=%4.2g$',...
        i,pcoords{:},ev),ltx{:});
    if ~pl1init
        pl1init=true;
        plseg=plot3(ax1,pcoords{:},'ko','MarkerSize',12,'MarkerFaceColor','r','HandleVisibility','off');
    else
        [plseg.XData,plseg.YData,plseg.ZData]=deal(pcoords{:});
    end
    pause(0.2);
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
    exportgraphics(figure(2),'ncsnic.pdf');
    saveas(figure(2),'nc_homsnic.fig');
end
