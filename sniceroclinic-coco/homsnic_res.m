function [data,res]=homsnic_res(~,data,y)
[iv,ip,funcs]=deal(data.iv,data.ip,data.funcs);
[    xeq1,      yeq1,      xeq2,      yeq2,      x0,      y0,      x1,      y1]=deal(...
y(iv.xeq1),y(iv.yeq1),y(iv.xeq2),y(iv.yeq2),y(iv.x0),y(iv.y0),y(iv.x1),y(iv.y1));
[    mu,      alpha,      beta,      gamma]=deal(...
y(iv.mu),y(iv.alpha),y(iv.beta),y(iv.gamma));
p([ip.mu,ip.alpha,ip.beta,ip.gamma])=[...
      mu;   alpha;   beta;   gamma];
[s1,s2]=deal(y(iv.s1),y(iv.s2));
p=p(:);
rsn=funcs.f([xeq1;yeq1],p);
[usn,usa,u0,u1]=deal([xeq1;yeq1],[xeq2;yeq2],[x0;y0],[x1;y1]);
Jsn=funcs.dfdx(usn,p);
rsn_d=det(Jsn);
rsa=funcs.f(usa,p);
% eigenvector for 0 in usn
vsn=evec(Jsn,u0-usn,@abs);
% eigenvector for stable ev in usa
vsa=evec(funcs.dfdx(usa,p),u1-usa,@(ev)ev);
% if svals are given check differences
if ~isnan(s1)
    rprojsn=usn+vsn*s1-u0;
    rprojsa=usa+vsa*s2-u1;
    res=[rsn;rsn_d;rsa;rprojsn;rprojsa];
    return
end
% else, suggest svals
res=[vsn'*(u0-usn);vsa'*(u1-usa)];
end
function v=evec(J,udir,srtfun)
[V,D]=eig(J);
ev=diag(D);
[~,iev]=sort(srtfun(ev));
v=V(:,iev(1));
v=v/norm(v)*sign(udir'*v);
end