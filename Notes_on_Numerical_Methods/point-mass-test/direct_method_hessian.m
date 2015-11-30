%
%  discretization of the constraints for the Direct Method Eqns. (29-34)
%
%  hessian objective(z) + sum lambda(k) * hessian constraint[k](z)
%
%
function jac = direct_method_hessian(z,sigma,lambda,auxdata)

  % step
  N     = auxdata.N ;
  nvars = auxdata.nvars ;
  h     = auxdata.h ;
  k2    = auxdata.k2 ;
  k3    = auxdata.k3 ;

  % non zeros elements of sparse jacobian
  nnz = 2*N+1 ;
  I   = zeros(1,nnz) ;
  J   = zeros(1,nnz) ;
  VAL = zeros(1,nnz) ;

  % z = [ x, v, u ], index offset for v(1) and u(1)
  start_v = N+1 ;

  i_start = 2*N+3 ;
  L = (h*k2*lambda(N+1:2*N) + k3*(lambda(i_start+1:i_start+N) + ...
                                  lambda(i_start+N+1:i_start+2*N)) )/2 ;

  % derivative v-v
  I(1)   = start_v+1 ;
  J(1)   = start_v+1 ;
  VAL(1) = L(1) ;
  nz     = 1 ;
  for k=2:N
    I(nz+1:nz+2)   = [ start_v+k,   start_v+k ] ;
    J(nz+1:nz+2)   = [ start_v+k-1, start_v+k ] ;
    VAL(nz+1:nz+2) = [ L(k-1), L(k)+L(k-1) ] ;
    nz             = nz + 2 ;
  end
  I(nz+1:nz+2)   = [ start_v+N+1, start_v+N+1 ] ;
  J(nz+1:nz+2)   = [ start_v+N,   start_v+N+1 ] ;
  VAL(nz+1:nz+2) = [ L(end), L(end) ] ;

  jac = sparse(I,J,VAL,nvars,nvars) ;
end
