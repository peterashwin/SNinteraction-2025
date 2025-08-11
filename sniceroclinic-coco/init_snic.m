function prob=init_snic(prob,cmd,snic_data,varargin)
default={'run','','lab',[],'isolargs',{},'table',[]};
opts=sco_set_options(default,varargin,'pass_on');
[iv,ip]=deal(snic_data.iv,snic_data.ip);
switch cmd
    case 'reload'
        reload=true;
        prob=coco_set(prob,'coll','MXCL','off');
        prob=ode_coll2coll(prob,'seg',opts.run,opts.lab);
    case {'init','start'}
        reload=false;
        funcs=snic_data.funcs;
        prob=ode_isol2coll(prob,'seg',funcs.f,funcs.dfdx,funcs.dfdp,opts.isolargs{:});
end
[uidxcoll,datacoll,u0]=coco_get_func_data(prob,'seg.coll','uidx','data','u0');
maps=datacoll.pr.coll_seg.maps;
idxold=[maps.p_idx;maps.x0_idx;maps.x1_idx];
prob=coco_add_pars(prob,'T',uidxcoll(maps.T_idx),'T');
prob=coco_add_pars(prob,'pars',uidxcoll(maps.p_idx),fieldnames(ip));
ivnames=fieldnames(iv);
newnames={'xeq','yeq','s1','s2'};%,'dist'};
if ~reload % find initial values for extra variables
    [ubc(:,1),ubc(:,2),p]=deal(u0(maps.x0_idx),u0(maps.x1_idx),u0(maps.p_idx));
    y0([iv.x0,   iv.y0,     iv.x1,      iv.y1])=[...
      ubc(1,1),ubc(2,1),ubc(1,end,1),ubc(2,end)];
    y0([iv.mu,   iv.alpha,   iv.beta,   iv.gamma])=[...
      p(ip.mu),p(ip.alpha),p(ip.beta),p(ip.gamma)];
    shvals=opts.table;
    y0( [iv.xeq,    iv.yeq,    iv.s1,    iv.s2])=...%,    iv.dist])=...
    [shvals.xeq,shvals.yeq,shvals.s1,shvals.s2];%,shvals.dist];
else
    ch=coco_read_solution('snic',opts.run,opts.lab,'chart');
    y0=ch.x;
end
y0=y0(:);
unew=y0(ismember(ivnames,newnames));
prob=coco_add_func(prob,'snic',@snic_res,snic_data,'zero','uidx',idxold,'u0',unew);
snic_idx=coco_get_func_data(prob,'snic','uidx');
prob=coco_add_pars(prob,'snic.pars',snic_idx([iv.xeq,iv.yeq,iv.s1,iv.s2,iv.x1,iv.y1]),...
    [newnames,{'x1','y1'}]);
prob=add_hetphasecond(prob,'seg.coll');
dist_idx=snic_idx([iv.xeq,iv.yeq,iv.x0,iv.y0,iv.x1,iv.y1]);
prob=coco_add_func(prob,'dist',@dist,[],'regular','dist','uidx',dist_idx);
end
function [data,val]=dist(~,data,y)
val=sqrt(sum((y(1:2)-y(3:4)).^2))+sqrt(sum((y(1:2)-y(5:6)).^2));
end