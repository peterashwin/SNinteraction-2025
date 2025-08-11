function [V0,W0,S,J]=snic_eigspace(df, xeq, par, x0)
%% extract right and left eigenvector for minimal-modulus eigenvalue 
J=df(xeq,par);
[U,S,V]=svd(J);
S=diag(S);
V0=V(:,end);
f0=x0-xeq;
V0=V0/(f0'*V0);
V0=V0/norm(V0);
W0=U(:,end);
W0=W0/(V0'*W0);
end