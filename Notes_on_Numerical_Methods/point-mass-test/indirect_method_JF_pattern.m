%% 
% map the indices with the corresponding index in the spase matrix
function jac = JfunPattern(auxdata)

  N     = auxdata.N ;
  nvars = auxdata.nvars ;

  % calcolo f(z)
  sx = 0 ;
  sv = sx+N+1 ;
  sl = sv+N+1 ;
  sm = sl+N+1 ;

  nnz = 16*N+4 ;
  I   = zeros(nnz,1) ;
  J   = zeros(nnz,1) ;

  eq = 0 ;
  nz = 0 ;
  % f(eq) = (x(k+1)-x(k))/h - vave(k) ;
  for k=1:N
    eq           = eq+1 ;
    I(nz+1:nz+4) = eq ;
    J(nz+1:nz+4) = [ sx+k+1, sx+k, sv+k+1,   sv+k ] ;
    nz           = nz+4 ;
  end
  % f(eq) = (v(k+1)-v(k))/h - u(k) + k0 + k1*vave(k) + k2*vave.^2
  for k=1:N
    eq           = eq+1 ;
    I(nz+1:nz+4) = eq ;
    J(nz+1:nz+4) = [ sm+k+1, sm+k, sv+k+1, sv+k ] ;
    nz           = nz+4 ;
  end
  % f(eq) = (lambda(k+1)-lambda(k))/h ;
  for k=1:N
    eq           = eq+1 ;
    I(nz+1:nz+2) = eq ;
    J(nz+1:nz+2) = [ sl+k+1, sl+k ] ;
    V(nz+1:nz+2) = [    1/h, -1/h ] ;
    nz           = nz+2 ;
  end
  % f(eq) = (mu(k+1)-mu(k))/h-1+lave-muave.*(k1+2*k2*vave)
  for k=1:N
    eq           = eq+1 ;
    I(nz+1:nz+6) = eq ;
    J(nz+1:nz+6) = [ sm+k+1, sm+k, sl+k+1, sl+k, sv+k+1, sv+k ] ;
    nz           = nz+6 ;
  end
  I(nz+1:nz+4) = [ eq+1, eq+2, eq+3, eq+4 ] ;
  J(nz+1:nz+4) = [ sx+1, sv+1, sv+N+1, sl+N+1 ] ;
  jac = sparse(I,J,ones(size(I)),nvars,nvars,nnz) ;
end