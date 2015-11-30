%% Compute NLP solution with IPOPT
%
% This function computes the numerical solution of the optimal control problem 
% as a direct transcription problem (NLP) and unsing IPOPT Optimizer.
% 
% See the paper for model details.
%
% Input parameters:
%  N                number of grid points
%  p_data           structure with problem data
%    p_data.g       acceleration of gravity g      ;
%    p_data.T_size  time horizon
%    p_data.k0      constant friction (normalized with mass)
%    p_data.k1      friction linearly dependent on speed (normalized with mass)
%    p_data.k2      drag friction (normalized with mass)
%    p_data.k3      down force (normalized with mass)
%  compute_accuracy compute Norm-2 w.r.t. exact solution: (true/false)
%                   Only for k2 = 0, k3 = 0
% plotting          plot solution: (true/false)
%
% Return:
%  norm2:          is empty if  compute_accuracy=false

function [z,elapsed,ok] = test_direct_method_with_ipopt(auxdata,use_hessian)

  [lb,ub,cl,cu] = direct_method_bound(auxdata) ; 

  options.lb = lb ; % Lower bound on the variables.
  options.ub = ub ; % Upper bound on the variables.

  % The constraint functions are bounded to zero
  options.cl = cl ;
  options.cu = cu ;

  % Set up the auxiliary data.
  options.auxdata = auxdata ;
  
  % Set the IPOPT options.
  options.ipopt.jac_d_constant   = 'no';
  options.ipopt.hessian_constant = 'yes';
  options.ipopt.mu_strategy      = 'adaptive';
  options.ipopt.max_iter         = 400;
  options.ipopt.tol              = 1e-10;
  
  % The callback functions.
  funcs.objective         = @direct_method_objective;
  funcs.gradient          = @direct_method_gradient;
  funcs.constraints       = @direct_method_constraints;
  funcs.jacobian          = @direct_method_constraints_jacobian;
  funcs.jacobianstructure = @direct_method_constraints_jacobian_pattern;
  if use_hessian
    funcs.hessian           = @direct_method_hessian;
    funcs.hessianstructure  = @direct_method_hessian_pattern;
  else
    options.ipopt.hessian_approximation      = 'limited-memory';
    options.ipopt.limited_memory_update_type = 'bfgs' ; % {bfgs}, sr1 = 6; % {6}
    options.ipopt.limited_memory_max_history = 50 ;
  end
  %options.ipopt.derivative_test = 'first-order';
  %options.ipopt.derivative_test = 'second-order';
  
  % Run IPOPT.
  z0 = direct_method_guess_solution(auxdata) ; 

  tic
  [z, info] = ipopt_auxdata(z0,funcs,options);
  elapsed = toc ;

  ok = info.iter ;

end
