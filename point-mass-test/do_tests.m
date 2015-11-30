%% Numerical Solution Method comparison for Optimal Control Problem
%
%  test different numerical solution methods for solving OCP problem.
%  Test compare the numerical methods vs computation time an on accuracy.
%  Accuracy is evaluated by solving the OCP problem over different grid 
%  and calculating the norm-2 of the error w.r.t. the analytical solution.

%% initialization
clc
close all
clear variables

addpath('../ipopt') ;
addpath('../scripts') ;

%  compute_accuracy compute Norm-2 w.r.t. exact solution: (true/false)
%                   Only for k2 = 0, k3 = 0
% plotting          plot solution: (true/false)
%
% Return:
%  norm2:          is empty if  compute_accuracy=false

g = 9.81;

% Data for non-linear test
p_data.g      = g      ; %  acceleration of gravity g  
p_data.T_size = 10     ; % time horizon
p_data.k0     = g*0.02 ; % constant friction (normalized with mass)
p_data.k1     = g*1e-5 ; % friction linearly dependent on speed (normalized with mass)
p_data.k2     = g*1e-4 ; % drag friction (normalized with mass)
p_data.k3     = g*2e-4 ; % down force (normalized with mass)

%grad_vec = [100,200,400,800,1000,1600,3200,6400,10000,12800];
grad_vec = [100,1000,10000];

%% Test numerical approach vs speed

comp_time = {} ;
kk        = 1 ;
for N=1000

  auxdata = direct_method_auxdata(N,p_data) ;
	
  if true
  	name           = sprintf('ipopt-no-hessian-%d',N) ;
    [z,elapsed,ok] = test_direct_method_with_ipopt(auxdata,false) ;
    norm2          = direct_method_save_solution(z,auxdata,name) ;
    comp_time{kk}  = sprintf('%25s elapsed = %14g, ||u-esatta|| = %g, converged = %d',name,elapsed,norm2,ok) ;
    kk = kk+1 ;
  end

  if true
  	name           = sprintf('ipopt-%d',N) ;
    [z,elapsed,ok] = test_direct_method_with_ipopt(auxdata,true);
    norm2          = direct_method_save_solution(z,auxdata,name) ;
    comp_time{kk}  = sprintf('%25s elapsed = %14g, ||u-esatta|| = %g, converged = %d',name,elapsed,norm2,ok) ;
    kk = kk+1 ;
  end

  if N <= 100
    name = sprintf('direct-matlab-%d',N) ;
    [z,elapsed,ok] = test_direct_method_with_matlab(auxdata);
    norm2          = direct_method_save_solution(z,auxdata,name) ;
    comp_time{kk}  = sprintf('%25s elapsed = %14g, ||u-esatta|| = %g, converged = %d',name,elapsed,norm2,ok) ;
    kk = kk+1 ;
  end
	
  % indirect

  auxdata = [] ;
  auxdata = indirect_method_auxdata(N,p_data) ;

  if N <= 10000
   	name           = sprintf('indirect-matlab-%d',N) ;
    [z,elapsed,ok] = test_indirect_method_with_matlab(auxdata);
    norm2          = indirect_method_save_solution(z,auxdata,name) ;
    comp_time{kk}  = sprintf('%25s elapsed = %14g, ||u-esatta|| = %g, converged = %d',name,elapsed,norm2,ok) ;
    kk = kk+1 ;
  end

  if true 
  	name           = sprintf('strscne-%d',N) ;
    [z,elapsed,ok] = test_indirect_method_with_strscne(auxdata);
    norm2          = indirect_method_save_solution(z,auxdata,name) ;
    comp_time{kk}  = sprintf('%25s elapsed = %14g, ||u-esatta|| = %g, converged = %d',name,elapsed,norm2,ok) ;
    kk = kk+1 ;
  end

  if true
  	name           = sprintf('tresnei-%d',N) ;
    [z,elapsed,ok] = test_indirect_method_with_tresnei(auxdata);
    norm2          = indirect_method_save_solution(z,auxdata,name) ;
    comp_time{kk}  = sprintf('%25s elapsed = %14g, ||u-esatta|| = %g, converged = %d',name,elapsed,norm2,ok) ;
    kk = kk+1 ;
  end

  if true
    name           = sprintf('affine_newton-%d',N) ;
    [z,elapsed,ok] = test_indirect_method_with_affine_newton(auxdata);
    norm2          = indirect_method_save_solution(z,auxdata,name) ;
    comp_time{kk}  = sprintf('%25s elapsed = %14g, ||u-esatta|| = %g, converged = %d',name,elapsed,norm2,ok) ;
    kk = kk+1 ;
  end
  
  indirect_method_plot_solution(z,auxdata,'affine_newton')
end

for k=1:length(comp_time)
  fprintf(1,'%s\n',comp_time{k}) ;    
end
