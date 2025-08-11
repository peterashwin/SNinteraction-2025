%------------------------
%% Define the system
%------------------------
%original equations are
%x_dot=x^2-lambda+2*alpha*y-gamma
%y_dot=y^2-lambda+2*beta*x+gamma

%curve of Hopf bifurcations/neutral saddles is
%lambda_tr0=gamma*(gamma+alpha^2-beta^2)/(alpha+beta)^2

%we set mu=lambda-lambda_tr0
params = {'mu', 'gamma', 'alpha', 'beta'};
ic=[params;num2cell(1:length(params))];
ip=struct(ic{:});
if exist('sys_hill', 'file')~=2
    syms x y
    syms(params)
    f=[x^2-(mu+gamma*(gamma+alpha^2-beta^2)/((alpha+beta)^2))+2*alpha*y-gamma;
        y^2-(mu+gamma*(gamma+alpha^2-beta^2)/((alpha+beta)^2))+2*beta*x+gamma];
    F = sco_sym2funcs(f, {[x; y], [mu; gamma; alpha; beta]}, ...
        {'x', 'p'}, 'maxorder', 3, 'filename', 'sys_hill');
else
    F=sco_gen(@sys_hill);
end
%assign the vector field and derivative to funcs
funcs = struct('f', F(''), 'dfdx', F('x'), 'dfdp', F('p'), ...
  'dfdx_dx', F('x*v'), 'dfdxdx_dx', F({'x','x*v'}), ...
  'dfdpdx_dx', F({'p','x*v'}));
thisbeta=-1.3;
