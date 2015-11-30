%%
%  discretization of dynamical system
%
function z = indirect_method_guess_solution(auxdata)
  % step
  T_size  = auxdata.T_size ;
  N       = auxdata.N ;
  g       = auxdata.g ;
  h       = T_size/N ;
  k0      = auxdata.k0 ;
  k1      = auxdata.k1 ;
  k2      = auxdata.k2 ;
  k3      = auxdata.k3 ;

  x      = zeros(N+1,1) ;
  v      = zeros(N+1,1) ;
  lambda = zeros(N+1,1) ;
  mu     = zeros(N+1,1) ;
  z      = [ x ; v ; lambda ; mu ] ;
end
