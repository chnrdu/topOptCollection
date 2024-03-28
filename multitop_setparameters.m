% as an example, the parameters are set as follows:
function [nx,ny,tol_out,tol_f,iter_max_in,iter_max_out,p,q,e,v,rf] = multitop_setparameters()
  nx = 100; ny = 50; tol_out = 0.001; tol_f = 0.05; iter_max_in = 2; iter_max_out = 400;
  p = 5; q = 3; e = [1000 400 200 100 1]'; v = [0.2 0.1 0.1 0.1 0.5]'; rf = 8;
end