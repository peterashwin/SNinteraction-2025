function [sol,u0,cseg]=coll_from_sol(collid,runid,lab)
[chart,data]=coco_read_solution(collid,runid,lab,'chart','data');
u0=chart.x;
cseg=data.coll_seg;
maps=cseg.maps;
par=u0(maps.p_idx);
t=cseg.mesh.tbp;
[degree,dim,ntst]=deal(cseg.int.NCOL,cseg.int.dim,cseg.maps.NTST);
xbp=reshape(u0(maps.xbp_idx),maps.xbp_shp);
sol=struct('profile',xbp,'parameter',par,'period',u0(maps.T_idx),...
    'mesh',reshape(t,degree+1,ntst),'degree',degree);