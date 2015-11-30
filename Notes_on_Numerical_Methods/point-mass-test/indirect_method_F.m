%
%
%
function f = indirect_method_F(z,auxdata)

  N  = auxdata.N ;
  h  = auxdata.h ;
  k0 = auxdata.k0 ;
  k1 = auxdata.k1 ;
  k2 = auxdata.k2 ;

  sx = 0 ;
  sv = sx+N+1 ;
  sl = sv+N+1 ;
  sm = sl+N+1 ;
  
  x      = z(sx+1:sv) ;
  v      = z(sv+1:sl) ;
  lambda = z(sl+1:sm) ;
  mu     = z(sm+1:sm+N+1) ;

  u     = indirect_method_u_eval(z,auxdata) ;

  vave  = (v(2:end)+v(1:end-1))/2 ;
  lave  = (lambda(2:end)+lambda(1:end-1))/2 ;
  muave = (mu(2:end)+mu(1:end-1))/2 ;

  f = [ (x(2:end)-x(1:end-1)) - h*vave ; ...
        (v(2:end)-v(1:end-1)) - h*(u - k0 - vave.*(k1+k2*vave)) ; ...
        (lambda(2:end)-lambda(1:end-1)) ; ...
        (mu(2:end)-mu(1:end-1))-h+h*lave-h*muave.*(k1+2*k2*vave) ; ...
        x(1) ; v(1) ; v(N+1) ; lambda(N+1) ] ;

end