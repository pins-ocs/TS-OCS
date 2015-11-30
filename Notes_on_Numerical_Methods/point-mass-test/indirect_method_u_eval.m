%% 
% map the indices with the corresponding index in the spase matrix
function u = indirect_method_u_eval(z,auxdata)

  g       = auxdata.g ;
  k3      = auxdata.k3 ;
  epsilon = auxdata.epsilon ;
  N       = auxdata.N ;

  % calcolo f(z)
  sx = 0 ;
  sv = sx+N+1 ;
  sl = sv+N+1 ;
  sm = sl+N+1 ;

  v  = z(sv+1:sl) ;
  mu = z(sm+1:sm+N+1) ;
  
  vave  = (v(1:end-1)+v(2:end))/2;
  muave = (mu(1:end-1)+mu(2:end))/2;

  tmp   = (2/pi)*(k3*vave.^2+g) ;
  u     = -tmp.*atan( (2/(pi*epsilon))*muave ) ;

end