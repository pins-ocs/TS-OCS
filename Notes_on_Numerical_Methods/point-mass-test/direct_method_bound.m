%%
%  discretization of dynamical system
%
function [lb,ub,cl,cu] = direct_method_bound(auxdata)
  N  = auxdata.N ;
  lb = [ -Inf*ones(2*(N+1),1) ; -Inf*ones(N,1) ] ;
  ub = [  Inf*ones(2*(N+1),1) ;  Inf*ones(N,1) ] ;
  cl = [ zeros(2*N+3,1) ; zeros(2*N,1) ] ; % 2*N+3 equality constraints
  cu = [ zeros(2*N+3,1) ; Inf*ones(2*N,1) ] ;
end
