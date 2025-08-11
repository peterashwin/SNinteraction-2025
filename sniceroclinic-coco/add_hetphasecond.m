function [prob, phasdata] = add_hetphasecond(prob, oid)
%% ad phase condition for heteroclinic
phasdata = init_data(prob, oid);
[fdata, uidx] = coco_get_func_data(prob, oid, 'data', 'uidx');
seg  = fdata.coll_seg;
maps = seg.maps;
prob = coco_add_func(prob, phasdata.fid, @FDF, phasdata, 'zero', ...
    'uidx', uidx(maps.xbp_idx), 'fdim', 1, ...
    'remesh', @remesh, 'F+DF');
prob = coco_add_slot(prob, phasdata.fid, @update, phasdata, 'update');
end

function phas = init_data(prob, oid)
%INIT_DATA   Initialize data for phase condition

[fdata, u0] = coco_get_func_data(prob,oid,'data','u0');

seg  = fdata.coll_seg;
maps = seg.maps;
mesh = seg.mesh;
phas = coco_func_data('protect');

phas.idx01  = [maps.x0_idx maps.x1_idx];
phas.xref   = u0(maps.xbp_idx);
phas.intfac = maps.Wp'*mesh.wts2*maps.W; % integral phase cond.
phas.pid=oid;
phas.fid    = coco_get_id(oid, 'phas');


end

function [phas, y, J] = FDF(prob, phas, u) %#ok<INUSL>
%FDF :phase condition

J0  =phas.xref'*phas.intfac;
xref=phas.xref;
[i0,i1]=deal(phas.idx01(:,1),phas.idx01(:,2));
[x1,x0,xr1,xr0]=deal(u(i1),u(i0),xref(i1),xref(i0));
y = J0*u+(x1/2-xr1)'*x1-(x0/2-xr0)'*x0;
J=J0;
J(i0)=J(i0)-(x0-xr0)';
J(i1)=J(i1)+(x1-xr1)';
end

function phas = update(prob, phas, cseg, varargin)

[fdata, uidx] = coco_get_func_data(prob, phas.pid, 'data', 'uidx');
seg  = fdata.coll_seg;
maps = seg.maps;

u         = cseg.src_chart.x(uidx);
phas.xref = u(maps.xbp_idx);
end

function [prob, status, xtr] = remesh(prob, phas, chart, old_u, old_V) %#ok<INUSD>
[fdata, uidx] = coco_get_func_data(prob, phas.pid, 'data', 'uidx');
seg  = fdata.coll_seg;
maps = seg.maps;
mesh = seg.mesh;
%u           = chart.x;
%phas.xref   = u(maps.xbp_idx);
phas.idx01  = [maps.x0_idx maps.x1_idx];
phas.intfac = maps.Wp'*mesh.wts2*maps.W; % integral phase cond.

xtr       = [];
prob      = coco_change_func(prob, phas, 'uidx', uidx(maps.xbp_idx), ...
  'fdim', 1);
status    = 'success';
end
