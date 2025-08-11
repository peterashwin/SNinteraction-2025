function  y = saddle_q(data, xbp, T0, T, par,funcs,ind,bd) %#ok<INUSD>
m    = data.coll_seg.int.NCOL;
xdim = data.coll_seg.int.dim;
tbp=data.coll_seg.mesh.tbp;
xbp=reshape(xbp,xdim,[]);
f = feval(funcs.f, xbp, repmat(par, [1 size(xbp, 2)])); % Evaluate vector field at basepoints
[fmin, idmin] = min(f(ind,:)); % Find basepoint closest to equilibrium
[fmax, idmax] = max(f(ind,:)); % Find basepoint closest to equilibrium
itrans=[min([idmin,idmax]),max([idmin,idmax])];
ttrans=tbp(itrans);
ip1=[1:itrans(1),itrans(2):size(xbp,2)];
ip2=itrans(1)+1:itrans(2)-1;
[fn1min,ieq1]=min(sum(f(:,ip1).*f(:,ip1),1));
ieq1=ip1(ieq1);
[fn2min,ieq2]=min(sum(f(:,ip2).*f(:,ip2),1));
ieq2=ip2(ieq2);
xeq1=xbp(:,ieq1);
t1=tbp(ieq1);
t2=tbp(ieq2);
xeq2=xbp(:,ieq2);
if fn1min>bd^2
    [xeq1out,det1,t1]=deal(NaN(size(xeq1)),NaN,NaN);
else
    xeq1out=xeq1;
    df=feval(funcs.dfdx, xeq1, par); % Evaluate Jacobian at equilibrium;
    det1=det(df);
end
if fn2min>bd^2
    [xeq2out,det2,t2]=deal(NaN(size(xeq2)),NaN,NaN);
else
    xeq2out=xeq2;
    df=feval(funcs.dfdx, xeq2, par); % Evaluate Jacobian at equilibrium;
    det2=det(df);
end
if xeq1(ind)<xeq2(ind)
    y=[ttrans;fn1min;xeq1out;t1;det1;fn2min;xeq2out;t2;det2];
else
    y=[ttrans([2,1]);fn2min;xeq2out;t2;det2;fn1min;xeq1out;t1;det1];
end
end