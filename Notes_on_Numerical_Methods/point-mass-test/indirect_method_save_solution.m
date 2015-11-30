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

function norm2 = indirect_method_save_solution(z,auxdata,name)

  % extract solution for plotting
  [x,v,u,lambda,mu,tt,xx,vv,uu] = indirect_method_extract_xvu(z,auxdata) ;
  
  t = (0:auxdata.h:auxdata.T_size).' ;
  N = auxdata.N ;

  % compute norm-2
  % filename for comparison with exact
  filename = sprintf('%s%d_e.txt',name,N);
  % check that k2 and k3 = 0. Exact solution is computeted under these
  % assumptions
  if auxdata.k2==0 && auxdata.k3==0
    auxdata
    %compute exact solution
    exact_sol = test_exact(N,auxdata,false);

    % compute norm-2 error
    norm2.x = norm(x-exact_sol.x)/sqrt(N+1) ;
    norm2.v = norm(v-exact_sol.v)/sqrt(N+1) ;
    norm2.u = norm(u-exact_sol.uc)/sqrt(N) ;
    %norm2.x = norm(x-exact_sol.x)/norm(exact_sol.x) ;
    %norm2.v = norm(v-exact_sol.v)/norm(exact_sol.v) ;
    %norm2.u = norm(u-exact_sol.uc)/norm(exact_sol.uc);
  else
    filename = sprintf('%s%d.txt',name,N);
    norm2 = -1;
  end
  
  % write out file
  ID = fopen(filename,'w') ;
  fprintf(ID,'t\tx\tv\tu\n') ;      
  for k=1:length(t)
    if k==1
      fprintf(ID,'%g\t%g\t%g\t%g\n',t(k),x(k),v(k),u(1)) ;
    else
      fprintf(ID,'%g\t%g\t%g\t%g\n',t(k),x(k),v(k),u(k-1)) ;
    end
  end
  fclose(ID) ;

end
