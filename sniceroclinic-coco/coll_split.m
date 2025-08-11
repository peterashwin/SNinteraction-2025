function [seg,eq,dev,S,V0,W0]=coll_split(sol,df)
%% shift large-period periodic solution s.t the approx equilbrium
% is at boundary, cut off mesh points to approximate distance dev
dsol=setfield(sol,'profile',coll_eva(sol,sol.mesh(:),'diff',1));
[dim,ncol1,ntst]=deal(size(sol.profile,1),size(sol.mesh,1),size(sol.mesh,2));
fnorm2=sum(dsol.profile.^2,1);
[~,imin]=min(fnorm2);
eq=sol.profile(:,imin);
tmin=sol.mesh(imin);
isplit=find(sol.mesh(1,:)<=tmin,1,'last');
x=reshape(sol.profile,[dim,ncol1,ntst]);
xshift=cat(3,x(:,:,isplit+1:end),x(:,:,1:isplit-1));
tshift=cat(2,sol.mesh(:,isplit+1:end),sol.mesh(end,end)+sol.mesh(:,1:isplit-1));
seg=setfield(sol,'profile',reshape(xshift,dim,ncol1*(ntst-1)));
seg.mesh=(tshift-tshift(1,1))/(tshift(end,end)-tshift(1,1));
seg.period=sol.period*(tshift(end,end)-tshift(1,1));
%% extract right and left eigenvector for minimal-modulus eigenvalue 
[V0,W0,S]=snic_eigspace(df, eq, sol.parameter, seg.profile(:,1));
dev=W0'*(seg.profile(:,[1,end])-eq(:,[1,1]));
dist=sqrt(sum((seg.profile(:,end)-eq).^2));
dev=[dev';dist];
end
