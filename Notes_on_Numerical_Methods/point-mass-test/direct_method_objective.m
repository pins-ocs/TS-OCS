%
%  discretization of the constraints for the Direct Method Eqns. (27)
%
%  sum (v[k+1]+v[k])/2
%
function f = direct_method_objective(z,auxdata)

	% split vector z into z, v, and u part Eqns. (28)
  [x,v,u] = direct_method_extract_xvu(z,auxdata) ;

  vave = (v(2:end)+v(1:end-1))/2 ;
  f    = -sum(vave) ;

end
