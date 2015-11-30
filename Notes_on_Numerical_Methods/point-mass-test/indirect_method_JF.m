%% 
% map the indices with the corresponding index in the spase matrix
function jac = Jfun(z,auxdata)

  N     = auxdata.N ;
  h     = auxdata.h ;
  k1    = auxdata.k1 ;
  k2    = auxdata.k2 ;
  nvars = auxdata.nvars ;

  % calcolo f(z)
  sx = 0 ;
  sv = sx+N+1 ;
  sl = sv+N+1 ;
  sm = sl+N+1 ;

  %x      = z(sx+1:sv) ;
  v      = z(sv+1:sl) ;
  %lambda = z(sl+1:sm) ;
  mu     = z(sm+1:sm+N+1) ;

  vave  = (v(2:end)+v(1:end-1))/2 ;
  muave = (mu(2:end)+mu(1:end-1))/2 ;

  DuDmu = indirect_method_u_Dmu(z,auxdata) ;
  DuDv  = indirect_method_u_Dv(z,auxdata) ;
  
  nnz = 16*N+4 ;
  I   = zeros(nnz,1) ;
  J   = zeros(nnz,1) ;
  V   = zeros(nnz,1) ;

  eq = 0 ;
  nz = 0 ;
  % f(eq) = (x(k+1)-x(k))/h - vave(k) ;
  for k=1:N
    eq           = eq+1 ;
    I(nz+1:nz+4) = eq ;
    J(nz+1:nz+4) = [ sx+k+1, sx+k, sv+k+1,   sv+k ] ;
    V(nz+1:nz+4) = [    1, -1, -(h/2), -(h/2) ] ;
    nz           = nz+4 ;
  end
  % f(eq) = (v(k+1)-v(k))/h - u(k) + k0 + vave.*(k1+k2*vave)
  for k=1:N
    eq           = eq+1 ;
    I(nz+1:nz+4) = eq ;
    J(nz+1:nz+4) = [ sm+k+1, sm+k, sv+k+1, sv+k ] ;
    V(nz+1:nz+4) = [ -h*DuDmu(k)/2, ...
                     -h*DuDmu(k)/2, ...
                      1 + h*k1/2 + h*k2*vave(k) - h*DuDv(k)/2, ...
                     -1 + h*k1/2 + h*k2*vave(k) - h*DuDv(k)/2 ] ;
    nz = nz+4 ;
  end
  % f(eq) = (lambda(k+1)-lambda(k))/h ;
  for k=1:N
    eq           = eq+1 ;
    I(nz+1:nz+2) = eq ;
    J(nz+1:nz+2) = [ sl+k+1, sl+k ] ;
    V(nz+1:nz+2) = [    1, -1 ] ;
    nz           = nz+2 ;
  end
  % f(eq) = (mu(k+1)-mu(k))/h-1+lave-muave.*(k1+2*k2*vave)
  for k=1:N
    eq           = eq+1 ;
    I(nz+1:nz+6) = eq ;
    J(nz+1:nz+6) = [ sm+k+1, sm+k, sl+k+1, sl+k, sv+k+1, sv+k ] ;
    V(nz+1:nz+6) = [  1-h*(k1/2+k2*vave(k)), ...
                     -1-h*(k1/2+k2*vave(k)), ...
                      h/2, ...
                      h/2, ...
                     -h*k2*muave(k), ...
                     -h*k2*muave(k) ] ;
    nz = nz+6 ;
  end
  %  x(1) ; v(1) ; v(N+1) ; lambda(N+1) ] ;
  I(nz+1:nz+4) = [ eq+1, eq+2, eq+3, eq+4 ] ;
  J(nz+1:nz+4) = [ sx+1, sv+1, sv+N+1, sl+N+1 ] ;
  V(nz+1:nz+4) = [ 1, 1, 1, 1 ] ;  
  
  jac = sparse(I,J,V,nvars,nvars,nnz) ;

end