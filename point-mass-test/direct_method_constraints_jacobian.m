%
%  discretization of the constraints for the Direct Method Eqns. (29-34)
%
%  (x[k+1]-x[k])/h-v[k+1/2] = 0
%  (v[k+1]-v[k])/h-u[k+1/2]-k0-k1*v[k+1/2]-k2*v[k+1/2]^2 = 0
%  x[0] = 0, v[0] = 0, v[N+1] = 0
%  g+k3*v[k+1/2]^2+u >= 0
%  g+k3*v[k+1/2]^2-u >= 0
%  Jacobian of cosntraints functions
%
%
function jac = direct_method_constraints_jacobian(z,auxdata)
  [~,v,~] = direct_method_extract_xvu(z,auxdata) ;

  % step
  %T_size = auxdata.T_size ;
  N      = auxdata.N ;
  g      = auxdata.g ;
  h      = auxdata.h ;
  nvars  = auxdata.nvars ;
  % k0     = auxdata.k0 ;
  k1     = auxdata.k1 ;
  k2     = auxdata.k2 ;
  k3     = auxdata.k3 ;

  %xave = (x(2:end)+x(1:end-1))/2 ;
  vave = (v(2:end)+v(1:end-1))/2 ;

  % non zeros elements of sparse jacobian
  nnz = 10*N+3 ;
  I   = zeros(1,nnz) ;
  J   = zeros(1,nnz) ;
  VAL = zeros(1,nnz) ;

  % z = [ x, v, u ], index offset for v(1) and u(1)
  start_v = N+1 ;
  start_u = start_v+N+1 ;

  %f = (x1 - x0)/h - vave ;
  neq = 0 ; % index of equation (row)
  nz  = 0 ;
  tmp = [ 1, -1, -(h/2), -(h/2) ] ; % row pattern
  for k=1:N
    neq = neq+1 ;
    I(nz+1:nz+4)   = neq ;
    J(nz+1:nz+4)   = [ k+1, k, start_v+k+1, start_v+k ] ;
    VAL(nz+1:nz+4) = tmp ;
    nz = nz + 4 ;
  end

  % (v1 - v0)/h - u + k0 + vave.*(k1+k2*vave) ; ...
  for k=1:N
    neq = neq+1 ;
    I(nz+1:nz+3)   = neq ;
    J(nz+1:nz+3)   = [ start_v+k+1, start_v+k, start_u+k ] ;
    VAL(nz+1:nz+3) = [ 1+h*(k1/2+k2*vave(k)), ...
                      -1+h*(k1/2+k2*vave(k)), ...
                      -h ] ;
    nz = nz + 3 ;
  end
	
  % boundary conditions
  %  x[0] = 0, v[0] = 0, v[N+1] = 0
  I(nz+1:nz+3)   = [ neq+1, neq+2,  neq+3 ] ;
  J(nz+1:nz+3)   = [ 1, start_v+1, start_v+N+1 ] ;
  VAL(nz+1:nz+3) = [ 1, 1,         1 ] ;
  neq = neq+3 ;
  nz  = nz+3 ;

  % g+k3*vave.^2+u
  tmp = k3*vave ;
  for k=1:N
    neq            = neq+1 ;
    I(nz+1:nz+3)   = neq ;
    J(nz+1:nz+3)   = [ start_v+k+1, start_v+k, start_u+k ] ;
    VAL(nz+1:nz+3) = [ tmp(k), tmp(k), 1 ] ;
    nz             = nz + 3 ;
  end
  % g+k3*vave.^2-u
  for k=1:N
    neq            = neq+1 ;
    I(nz+1:nz+3)   = neq ;
    J(nz+1:nz+3)   = [ start_v+k+1, start_v+k, start_u+k ] ;
    VAL(nz+1:nz+3) = [ tmp(k), tmp(k), -1 ] ;
    nz             = nz + 3 ;
  end

  jac = sparse(I,J,VAL,neq,nvars) ;
end
