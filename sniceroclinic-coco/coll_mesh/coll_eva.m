function out=coll_eva(sol,x,varargin)
dim=size(sol.profile,1);
default={'kron',dim,'output','profile','dims',1:dim};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
if strcmp(options.output,'profile')
    options.kron=length(options.dims);
end
submesh=sol.mesh(:,1);
submesh=2*(submesh-submesh(1))/(submesh(end)-submesh(1))-1;
msh=[sol.mesh(1,:),sol.mesh(end,end)];
J=coll_mesh_mat(msh,submesh,x,pass_on{:},'kron',options.kron);
if strcmp(options.output,'profile')
    out=J*reshape(sol.profile(options.dims,:),[],1);
    out=reshape(out,[length(options.dims),size(x)]);
else
    out=J;
end
end
