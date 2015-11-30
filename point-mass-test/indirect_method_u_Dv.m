%% 
% map the indices with the corresponding index in the spase matrix
function DuDv = indirect_method_u_Dv(z,auxdata)

  g       = auxdata.g ;
  k3      = auxdata.k3 ;
  epsilon = auxdata.epsilon ;
  N       = auxdata.N ;
  
  % calcolo f(z)
  sx = 0 ;
  sv = sx+N+1 ;
  sl = sv+N+1 ;
  sm = sl+N+1 ;
  su = sm+N+1 ;

  v  = z(sv+1:sl) ;
  mu = z(sm+1:su) ;
  
  vave  = (v(1:end-1)+v(2:end))/2;
  muave = (mu(1:end-1)+mu(2:end))/2;
  DuDv  = -(4*k3/pi)*(atan((2/(pi*epsilon))*muave).*vave) ;

  %for k=1:N
  %  vave  = (v(k)+v(k+1))/2 ;
  %  muave = (mu(k)+mu(k+1))/2 ;
  %  DuDv(k) = -4*atan(2*muave/(pi*epsilon))*a3*vave/pi ;
  %end

end