%
%  discretization of the constraints for the Direct Method Eqns. (27)
%
%  (x[k+1]-x[k])/h-v[k+1/2] = 0
%  (v[k+1]-v[k])/h-u[k+1/2]-k0-k1*v[k+1/2]-k2*v[k+1/2]^2 = 0
%  x[0] = 0, v[0] = 0, v[N+1] = 0
%  -1 <= u/(g+k3*v[k+1/2]^2) <= 1
%
function f = direct_method_constraints(z,auxdata)

  % split vector z into z, v, and u part Eqns. (28)
  [x,v,u] = direct_method_extract_xvu(z,auxdata) ;

  % extract aux parameters
  g  = auxdata.g ;
  h  = auxdata.h ;
  k0 = auxdata.k0 ;
  k1 = auxdata.k1 ;
  k2 = auxdata.k2 ;
  k3 = auxdata.k3 ;

  vave   = (v(2:end)+v(1:end-1))/2 ;
  deltax = x(2:end)-x(1:end-1) ;
  deltav = v(2:end)-v(1:end-1) ;

  f = [ deltax - h*vave ; ...
        deltav - h*(u - k0 - vave.*(k1+k2*vave)) ; ...
        x(1) ; v(1) ; v(end) ; ...
        (g+k3*vave.^2)+u ; ...
        (g+k3*vave.^2)-u ] ;
end
