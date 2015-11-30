%%
%  Setup data structure for computation
%
function auxdata = direct_method_auxdata(N,p_data)
  
  nvars          = 3*N+2  ;
  nconts         = 4*N+3  ; % x, v, u min, u max, BC
  auxdata.N      = N      ;
  auxdata.nvars  = nvars  ;
  auxdata.nconts = nconts ;
  auxdata.g      = p_data.g ;
  auxdata.T_size = p_data.T_size ; % final time
  auxdata.h      = auxdata.T_size/auxdata.N ;
  auxdata.k0     = p_data.k0; 
  auxdata.k1     = p_data.k1; 
  auxdata.k2     = p_data.k2;
  auxdata.k3     = p_data.k3;
  
end
