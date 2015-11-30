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

function direct_method_plot_solution(z,auxdata,name)

  % extract solution for plotting
  [x,v,u,tt,xx,vv,uu] = direct_method_extract_xvu(z,auxdata) ;
  
  t = (0:auxdata.h:auxdata.T_size).' ;
  
  % visualize solution
  subplot(3,1,1) ;
  plot(tt,uu) ;
  title('control') ;
    
  subplot(3,1,2) ;
  plot(tt,xx) ;
  title('space') ;

  subplot(3,1,3) ;
  plot(tt,vv) ;
  title('velocity') ;

end
