function animate_timeprofiles(iv,ip,funcs,bdrun,snbd,hombd,varargin)
default={'pause',{0.2},'figure',2,'h0',1e-1,'labels',[]};
opts=dde_set_options(default,varargin,'pass_on');
bdsnic=coco_bd_table(bdrun);
labs=bdsnic{:,'LAB'};
if iscell(labs)
    labrows=cellfun(@(s)~isempty(s),labs);
    labs=cell2mat(bdsnic{labrows,'LAB'});
end
figure(opts.figure);clf;
tl=tiledlayout(1,2);nexttile;ax1=gca;
lw={'linewidth',2};
ltx={'Interpreter','latex','FontSize',20};
txt={'FontSize',20,'FontName','Courier','FontWeight','bold'};
hold(ax1,'on');
clr=lines();
plot(ax1,hombd.mu,hombd.gamma,'.-',snbd.mu,snbd.gamma,'.-');
set(ax1, txt{:},lw{:},'box','on')
xlim(ax1,[-0.12,0.03]);
ylim(ax1,[3.54,3.66]);
pl1init=false;
nexttile;ax2=gca;
h0=opts.h0;
if isempty(opts.labels)
    irg=1:length(labs);
end
for i=irg
    lab=labs(i);
    try
        [cs,ds]=coll_read_solution('seg',bdrun,lab,'chart','data');
    catch
        continue
    end
    t=cs.tbp;
    u=cs.xbp;
    ch=coco_read_solution('snic',bdrun,lab,'chart');
    y0=ch.x;
    plot(ax2,u(:,1),u(:,2),'o-',lw{:});
    hold(ax2,'on');
    usn=[y0(iv.xeq);y0(iv.yeq)];
    plot(ax2,usn(1),usn(2),'ko','MarkerSize',10,'MarkerFaceColor','y');
    vsn=snic_eigspace(funcs.dfdx,usn,cs.p,y0([iv.x0,iv.y0]));
    plot(ax2,usn(1)+h0*vsn(1)*[-1,1],usn(2)+h0*vsn(2)*[-1,1],'k-',lw{:});
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
        i,pcoords{:},ch.x(iv.s1),ch.x(iv.s2)),ltx{:});
    if ~pl1init
        pl1init=true;
        plseg=plot(ax1,pcoords{1:2},'ko','MarkerSize',12,'MarkerFaceColor','r','HandleVisibility','off');
    else
        [plseg.XData,plseg.YData]=deal(pcoords{1:2});
    end
    pause(opts.pause{:});
    drawnow
end
end