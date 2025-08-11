%% test dde_coll_mesh_mat
clear
deg=35;   % degree of interpolation polynomials on each subinterval
nint=400; % number of subintervals
rng(1);   % random number seed for repeatability
T0=-pi/2;
T=2*pi;   % interval [T0,T0+T] on which collocation is done
coarse_msh=[0,cumsum(rand(1,nint))]; % collocation mesh with random subintervals
coarse_msh=coarse_msh/coarse_msh(end)*T+T0; %adjusted to [T0,T0+T] 
submesh=dde_coll_set_grid('storage',deg); % Chebyshev nodes, 2nd kind, by default
%% create equidistant submesh (by default)
m_eq=dde_coll_meshfill(coarse_msh,deg);
%% Chebyshev nodes 2nd kind
m_ch=dde_coll_meshfill(coarse_msh,deg,'grid','cheb');
%% Gauss-Legendre nodes for collocation
m_gl=dde_coll_meshfill(coarse_msh,deg,'grid','gauss');
%% interpolation mapping equidistant mesh onto Cheb mesh
J_ch_eq=coll_mesh_mat(coarse_msh,linspace(-1,1,deg+1),m_ch);
%% interpolation mapping Cheb mesh onto Legendre mesh
J_gl_ch=coll_mesh_mat(coarse_msh,submesh,m_gl);
%% interpolation with shift and wrap, mapping Cheb mesh onto Legendre mesh
tau=1;
wrap=@(m)mod(m-T0,T)+T0;
J_gltau_ch=coll_mesh_mat(coarse_msh,submesh,wrap(m_gl-tau));
%% Differentiation with shift and wrap, mapping Cheb mesh onto Legendre mesh
tau=1;
J_Dgltau_ch=coll_mesh_mat(coarse_msh,submesh,wrap(m_gl-tau),'diff',1);
%% Same for function profiles of higher dimension (kron)
J2_Dgltau_ch=coll_mesh_mat(coarse_msh,submesh,wrap(m_gl-tau),'diff',1,'kron',2);
%% tests
vec=@(x)x(:);
n=@(x)norm(x,'inf');
err_ch_eq=n(J_ch_eq*vec(sin(m_eq))-vec(sin(m_ch))) %#ok<*NOPTS>
err_gl_ch=n(J_gl_ch*vec(sin(m_ch))-vec(sin(m_gl)))
err_gltau_ch=n(J_gltau_ch*vec(sin(m_ch))-vec(sin(m_gl-tau)))
err_Dgltau_ch=n(J_Dgltau_ch*vec(sin(m_ch))-vec(cos(m_gl-tau)))
err2_Dgltau_ch=n(...
    J2_Dgltau_ch*vec([sin(m_ch);cos(m_ch)])-...
    vec([cos(m_gl-tau);-sin(m_gl-tau)]))
%% structure of interpolation matrix
spy(J2_Dgltau_ch)
title('structure of $W''(-\tau)$','Interpreter','LaTeX')