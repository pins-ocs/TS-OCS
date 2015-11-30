%
%  discretization of the constraints for the Direct Method Eqns. (27)
%
%  sum (v[k+1]+v[k])/2
%
function grad = direct_method_gradient(z,auxdata)

	% split vector z into z, v, and u part Eqns. (28)
  [x,v,u] = direct_method_extract_xvu(z,auxdata) ;
  N      = auxdata.N ;
  nvars  = auxdata.nvars ;

  % z = [ x, v, u ], index offset for v(1) and u(1)
  start_v = N+1 ;
  start_u = start_v+N+1 ;

  grad                = zeros(nvars,1) ;
  grad(start_v+1)     = -1/2 ;
  grad(start_v+(2:N)) = -1   ;
  grad(start_v+N+1)   = -1/2 ;

end
