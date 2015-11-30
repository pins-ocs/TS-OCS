%%
%  discretization of dynamical system
%
function [x,v,u,varargout] = direct_method_extract_xvu(z,auxdata)
  % map x, v, lambda, mu, u from vector z
  N      = auxdata.N ;
  T_size = auxdata.T_size ;
  h      = T_size/N ;
	
	start_v = N+1 ;
	start_u = start_v+N+1 ;

  x = z(1:N+1) ;
  v = z(start_v+(1:N+1)) ;
  u = z(start_u+(1:N)) ;

  if nargout == 7
    T_size = auxdata.T_size ;
    N      = auxdata.N ;
    h      = T_size/N ;
    tt = zeros(2*N,1) ;
    xx = zeros(2*N,1) ;
    vv = zeros(2*N,1) ;
    uu = zeros(2*N,1) ;
    for k=1:N
      tt(2*k-1) = (k-1)*h ;
      tt(2*k)   = k*h ;
      xx(2*k-1) = x(k) ;
      xx(2*k)   = x(k+1) ;
      vv(2*k-1) = v(k) ;
      vv(2*k)   = v(k+1) ;
      uu(2*k-1) = u(k) ;
      uu(2*k)   = u(k) ;
    end
    varargout{1} = tt ;
    varargout{2} = xx ;
    varargout{3} = vv ;
    varargout{4} = uu ;
  end
  % u = 0.5*([0 ; u] + [ u ; 0 ]) ;
end