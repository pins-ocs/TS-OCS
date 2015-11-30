%% 
% map the indices with the corresponding index in the spase matrix
function varargout = indirect_method_f_model(z)
  global auxdata ;
  varargout{1} = indirect_method_F(z,auxdata) ;
  if nargout > 1
    varargout{2} = indirect_method_JF(z,auxdata) ;
  end
end