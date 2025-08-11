function [data,res]=snic_res(~,data,y)
[iv,ip,funcs]=deal(data.iv,data.ip,data.funcs);
[    xeq,      yeq,        x0,      y0,      x1,      y1]=deal(...
y(iv.xeq),y(iv.yeq),y(iv.x0),y(iv.y0),y(iv.x1),y(iv.y1));
[    mu,      alpha,      beta,      gamma]=deal(...
y(iv.mu),y(iv.alpha),y(iv.beta),y(iv.gamma));
p([ip.mu,ip.alpha,ip.beta,ip.gamma])=[...
      mu;   alpha;   beta;   gamma];
p=p(:);
[s1,s2]=deal(y(iv.s1),y(iv.s2));%,y(iv.dist));
[usn,ubc(:,1),ubc(:,2)]=deal([xeq;yeq],[x0;y0],[x1;y1]);
% right and left eigenvector for 0 in usn
[vsn,wsn,~,Jsn]=snic_eigspace(funcs.dfdx,usn,p,ubc(:,1));
rsn=funcs.f([xeq;yeq],p);
rsn_d=det(Jsn);
% if svals are given check differences
if ~isnan(s1)
    rprojsn=usn+vsn*s1-ubc(:,1);
    lprojsn=wsn'*(ubc(:,2)-usn)-s2;
    %resdist=sum((ubc(:,2)-usn).^2)-dist^2;
    res=[rsn;rsn_d;rprojsn;lprojsn];%;resdist];
    return
end
% else, suggest svals
res=[wsn'*(ubc(:,1)-usn);wsn'*(usn-ubc(:,2));sqrt(sum(ubc(:,2)-usn).^2)];
end
