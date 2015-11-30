%%
%  discretization of dynamical system
%
function jac = direct_method_constraints_jacobian_pattern(auxdata)

  % step
  T_size = auxdata.T_size ;
  N      = auxdata.N ;
  nvars  = auxdata.nvars ;

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
  for k=1:N
    neq = neq+1 ;
    I(nz+1:nz+4) = neq ;
    J(nz+1:nz+4) = [ k+1, k, start_v+k+1, start_v+k ] ;
    nz = nz + 4 ;
  end

  % (v1 - v0)/h - u + k0 + vave.*(k1+k2*vave) ; ...
  for k=1:N
    neq = neq+1 ;
    I(nz+1:nz+3) = neq ;
    J(nz+1:nz+3) = [ start_v+k+1, start_v+k, start_u+k ] ;
    nz = nz + 3 ;
  end
	
	% boundary conditions
	%  x[0] = 0, v[0] = 0, v[N+1] = 0
  I(nz+1:nz+3)   = [ neq+1, neq+2,  neq+3 ] ;
  J(nz+1:nz+3)   = [ 1, start_v+1, start_v+N+1 ] ;
  neq = neq+3 ;
  nz  = nz+3 ;

  % g+k3*vave.^2+u
  for k=1:N
    neq = neq+1 ;
    I(nz+1:nz+3) = neq ;
    J(nz+1:nz+3) = [ start_v+k+1, start_v+k, start_u+k ] ;
    nz = nz + 3 ;
  end
  % g+k3*vave.^2-1
  for k=1:N
    neq = neq+1 ;
    I(nz+1:nz+3) = neq ;
    J(nz+1:nz+3) = [ start_v+k+1, start_v+k, start_u+k ] ;
    nz = nz + 3 ;
  end

  jac = sparse(I,J,ones(1,nz),neq,nvars) ;

end
