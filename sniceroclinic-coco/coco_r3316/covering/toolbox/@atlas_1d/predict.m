function [prob, atlas, cseg, correct] = predict(atlas, prob, cseg) %#ok<INUSD>
chart = atlas.chart_list{1};
if ~atlas.first && atlas.cont.NAdapt*chart.pt > 0 ...
    && mod(chart.pt,atlas.cont.NAdapt) == 0
  atlas.first = true;
  % remesh base chart of new curve segment
  
  x0        = chart.x;
  V0        = [ chart.t chart.TS ];
  [prob, chart, x0, V0] = coco_remesh(prob, chart, x0, V0, atlas.cont.RMMX);
  nv        = repmat(sqrt(sum(V0.^2,1)), size(V0,1), 1);
  chart.x   = x0;
  chart.t   = V0(:,1)./nv(:,1);
  chart.TS  = V0(:,2:end)./nv(:,2:end);
  
  % update pseudo-arclength projection condition
  prcond.x  = chart.x;
  prcond.TS = chart.TS;
  prcond.s  = chart.s;
  prcond.h  = 0;
  x1        = chart.x;
  % construct new curve segment including potential problem update
  [prob, cseg] = CurveSegment.create(prob, chart, prcond, x1);
  
  % update NV for non-square problems
  [prob, chart, J] = prob.efunc.DFDX(prob, chart, x1); %#ok<ASGLU>
  [m, n] = size(J);
  if m>n
    [L, U, P]  = lu(J); %#ok<ASGLU>
    Y  = L(1:n, 1:n)' \ L(n+1:end, 1:n)';
    B(:,1:m-n) = P'*[Y; -speye(m-n)];
    cseg.curr_chart.NV = orth(full(B));
  end
else
  atlas.first = false;
  
  % update pseudo-arclength projection condition
  prcond.x  = chart.x;
  prcond.TS = chart.TS;
  prcond.s  = chart.s;
  prcond.h  = chart.R/norm(chart.TS, atlas.cont.norm);
  x1        = chart.x + chart.TS*(prcond.h*chart.s);
  % construct new curve segment
  [prob, cseg] = CurveSegment.create(prob, chart, prcond, x1);
end
correct = true;

% bug: to re-implement option log_angle functionality get old code with:
% svn cat -r2839 atlas_1d.m | tail -n +313 | head -62
end
