%% Evaluate (derivative of) piecewise collocation polynomial
% assuming all requested points are inside the collocation interval
% mesh (assertion is placed)
%
% $Id$
%%
function J=coll_mesh_mat(msh,submesh,x,varargin)
%% INPUT:
%
% * msh: 1x(N+1) sequence of strictly monotone increasing times
% * submesh: m base points, scaled to [-1,1], for polynomial
% * x (1 x nx): point(s) where to evaluate
%
%% Optional
%
% * loc of size(1,nx)): same size as x: if x(k)==msh(i) for some i in 2:N,
% then x fits to the interval msh([i-1,i]) if loc(k)==-1,
% otherwise, if loc(k)==+1, x fits to the interpolation in msh([i,i+1])
% For all other x loc is irrelevant. Default for loc is that for any i such
% that x(i)==x(i+1), loc(i) is -1, for all others it is 1.
% * kron (1): if kron(J,speye(n)) should be done at the end
% * sparse (true): matrix returned should be sparse
%
%% OUTPUT: 
%   J: sparse Jacobian:  P(y(:))(:)=J*y(:) for the piecewise collocation
%   polynomial on fullmesh(:)
%% default for loc
x=x(:)';
assert(all(x>=msh(1))&&all(x<=msh(end)));
nx=length(x);
locdefault=ones(1,numel(x));
locdefault([x(1:end-1)==x(2:end),false])=-1;
%% options
default={'kron',1,'diff',0,'sparse',true,'loc',locdefault};
options=dde_set_options(default,varargin,'pass_on');
%% fill in complete mesh
msh=msh(:)';
nt=length(msh);
nint=nt-1;
submesh=submesh(:);
deg1=length(submesh);
degree=deg1-1;
dt=diff(msh);
%% locate in which subinterval each x is (ix)
ix=floor(interp1(msh,1:nt,x,'linear'));
%% adjust to lower interval if requested or required at final boundary
l_adj= (x==msh(ix) & options.loc(:)'==-1 & ix>1) | ix==nt;
ix(l_adj)=ix(l_adj)-1;
%% evaluate Lagrange polynomials on all interpolation times
xscal=2*(x-msh(ix))./dt(ix)-1;
%% evaluate barycentric weights
% for flexibility, (not assuming any particular interpolation grid)
[w,Dmat]=barywt(submesh,options.diff);
%% calculate interpolation matrix
% using barycentric interpolation
ti_m=(ix-1)*(degree+1)+1;
og=ones(deg1,1);
ox=ones(nx,1);
denom=xscal(og,:)-submesh(:,ox);
jac_ind(:,:,2)=ti_m(og,:)+repmat((0:degree)',1,nx);
jac_ind(:,:,1)=repmat(1:nx,deg1,1);
fac=w(:,ox)./denom;
jac_vals=zeros(deg1,nx);
jac_vals(~denom(:))=1;
denomfin=all(denom~=0,1);
jac_vals(:,denomfin)=fac(:,denomfin)./repmat(sum(fac(:,denomfin),1),deg1,1);
%% Differentiate if requested
if options.diff>0
    div=dt(ix).^options.diff;
    jac_vals=(Dmat'*jac_vals)./repmat(div,deg1,1);
end
%% assemble Jacobian
jac_ind=reshape(jac_ind,[],2);
jac_vals=jac_vals(:);
J=sparse(jac_ind(:,1),jac_ind(:,2),jac_vals,nx,deg1*nint);
if options.kron>1
    J=kron(J,speye(options.kron));
end
if ~options.sparse
    J=full(J);
end
end
%% construct barycentric weights for mesh
function [w,Dout]=barywt(submesh,difforder)
npi=length(submesh);
base_h=submesh';
base_v=submesh;
xdiff=base_h(ones(npi,1),:)-base_v(:,ones(npi,1));
w=1./prod(eye(npi)-xdiff,2);
Dout=eye(npi);
if difforder>0
    wrep=w(:,ones(npi,1));
    denom=(xdiff'-xdiff);
    denom(1:npi+1:end)=Inf;
    D=wrep'./wrep./denom;
    Drsum=sum(D,2);
    D(1:npi+1:end)=-Drsum;
    for i=1:difforder
        Dout=4*D*Dout;
    end
end
end