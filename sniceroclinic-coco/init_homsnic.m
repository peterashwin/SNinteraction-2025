function prob=init_homsnic(prob,cmd,homsnic_data,varargin)
default={'run','','lab',[],'isolargs',{},'table',[]};
opts=sco_set_options(default,varargin,'pass_on');
[iv,ip]=deal(homsnic_data.iv,homsnic_data.ip);
switch cmd
    case 'reload'
        reload=true;
        prob=coco_set(prob,'coll','MXCL','off');
        prob=ode_coll2coll(prob,'seg',opts.run,opts.lab);
    case {'init','start'}
        reload=false;
        funcs=homsnic_data.funcs;
        prob=ode_isol2coll(prob,'seg',funcs.f,funcs.dfdx,funcs.dfdp,opts.isolargs{:});
end
[uidxcoll,datacoll,u0]=coco_get_func_data(prob,'seg.coll','uidx','data','u0');
maps=datacoll.pr.coll_seg.maps;
idxold=[maps.p_idx;maps.x0_idx;maps.x1_idx];
prob=coco_add_pars(prob,'T',uidxcoll(maps.T_idx),'T');
ivnames=fieldnames(iv);
newnames={'xeq1','yeq1','xeq2','yeq2','s1','s2'};
if ~reload % find initial values for extra variables
    [ubc(:,1),ubc(:,2),p]=deal(u0(maps.x0_idx),u0(maps.x1_idx),u0(maps.p_idx));
    y0([iv.x0,   iv.y0,     iv.x1,      iv.y1])=[...
      ubc(1,1),ubc(2,1),ubc(1,end,1),ubc(2,end)];
    y0([iv.mu,   iv.alpha,   iv.beta,   iv.gamma])=[...
      p(ip.mu),p(ip.alpha),p(ip.beta),p(ip.gamma)];
    shvals=opts.table;
    y0([iv.xeq1,    iv.yeq1,    iv.xeq2,    iv.yeq2])=[...
    shvals.xeq1,shvals.yeq1,shvals.xeq2,shvals.yeq2];
    y0([iv.s1,iv.s2])=NaN;
    y0=y0(:);
    [~,y0([iv.s1,iv.s2])]=homsnic_res(prob,homsnic_data,y0);
    for i=length(newnames):-1:1
        unew(i)=y0([iv.(newnames{i})]);
    end
    unew=unew(:);
else
    ch=coco_read_solution('homsnic',opts.run,opts.lab,'chart');
    unew=ch.x(ismember(ivnames,newnames));
end
prob=coco_add_func(prob,'homsnic',@homsnic_res,homsnic_data,'zero','uidx',idxold,'u0',unew);
homsnic_idx=coco_get_func_data(prob,'homsnic','uidx');
prob=coco_add_pars(prob,'homsnic.pars',homsnic_idx([iv.xeq1,iv.yeq1,iv.xeq2,iv.yeq2,iv.s1,iv.s2]),newnames);
end