function [z,elapsed,ok] = test_indirect_method_with_affine_newton(auxdata_in)

  global auxdata ;
  addpath('../scripts') ;
  
  auxdata = auxdata_in ;
 
  % compute guess solution
  z0 = indirect_method_guess_solution( auxdata ) ;
  ls = 1 ;

  tic
  tol  = 1e-10 ;
  opts = {'linesearch', ls, 'tol', 1E-10, 'maxiter', 100} ;
  auxdata.epsilon = 1e-1 ;
  [z1,f1,nit,ierr] = NewtonNonlinear( @indirect_method_f_model, z0, opts{:} ) ;
  ok = nit ;
  if ierr == 0
    auxdata.epsilon = 1e-3 ;
    [z2,f1,nit,ierr] = NewtonNonlinear( @indirect_method_f_model, z1, opts{:} ) ;
    ok = ok + nit ;
    auxdata.epsilon = auxdata_in.epsilon ;
    if ierr == 0
      [z,f1,nit,ierr] = NewtonNonlinear( @indirect_method_f_model, z2, opts{:} ) ;
      ok = ok + nit ;
    end
  end
  elapsed = toc ;
  
  if ierr ~= 0
    ok = -1 ;
  end

end
