function [G,diff_mat]=gmmjacob_num(f,coef,x,param)

diff_mat=numeric_jacob(@(param)f(x,param), param);
G=coef'*diff_mat;

end