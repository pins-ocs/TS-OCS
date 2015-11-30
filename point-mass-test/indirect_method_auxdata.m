%%
%  Setup data structure for computation
%
function auxdata = indirect_method_auxdata(N,p_data)
  
  auxdata.N       = N ;
  auxdata.nvars   = 4*N+4  ;
  auxdata.g       = p_data.g ;
  auxdata.T_size  = p_data.T_size ; % final time
  auxdata.h       = auxdata.T_size/auxdata.N ;
  auxdata.k0      = p_data.k0 ;
  auxdata.k1      = p_data.k1 ; 
  auxdata.k2      = p_data.k2 ;
  auxdata.k3      = p_data.k3 ;
	auxdata.epsilon = 1e-6 ;

end
