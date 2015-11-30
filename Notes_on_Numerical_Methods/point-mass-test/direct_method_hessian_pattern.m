%
%  discretization of the constraints for the Direct Method Eqns. (29-34)
%
%  hessian objective(z) + sum lambda(k) * hessian constraint[k](z)
%
%  The hessian of objective is zero
%
%  The hessian of constaints: (x[k+1]-x[k])/h-v[k+1/2] is 0
%  The hessian of constaints: (v[k+1]-v[k])/h-u[k+1/2]-k0-k1*v[k+1/2]-k2*v[k+1/2]^2 = 0
%  is the second derivative of -k2*v[k+1/2]^2
%  The hessian of constaints:  x[0] = 0, v[0] = 0, v[N+1] = 0 is zero
%  The hessian of constaints:  g+k3*v[k+1/2]^2+u
%  The hessian of constaints:  g+k3*v[k+1/2]^2-u
%
function jac = direct_method_hessian_pattern(auxdata)

  %T_size = auxdata.T_size ;
  N     = auxdata.N ;
  nvars = auxdata.nvars ;

  % non zeros elements of sparse jacobian
  nnz = 2*N+1 ;
  I   = zeros(1,nnz) ;
  J   = zeros(1,nnz) ;

  % z = [ x, v, u ], index offset for v(1) and u(1)
  start_v = N+1 ;

  % lambda.*k2*vave.^2 ; ...
  I(1) = start_v+1 ;
  J(1) = start_v+1 ;
  nz   = 1 ;
  for k=2:N
    I(nz+1:nz+2) = [ start_v+k,   start_v+k ] ;
    J(nz+1:nz+2) = [ start_v+k-1, start_v+k ] ;
    nz           = nz + 2 ;
  end
  I(nz+1:nz+2) = [ start_v+N+1, start_v+N+1 ] ;
  J(nz+1:nz+2) = [ start_v+N,   start_v+N+1 ] ;
  nz           = nz + 2 ;
  jac = sparse(I,J,ones(nz,1),nvars,nvars) ;
end
