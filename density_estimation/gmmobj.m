function [val,grad] = gmmobj(param,x,coef,moment,weight,fun_pdf)

diff=(fun_pdf(x,param)'*coef)'-moment;
% Braess
% diff=((ones(size(x))*1/param.*(x<=param))'*coef)'-moment;
val=diff'*weight*diff;

grad=2*diff'*weight*gmmjacob_num(fun_pdf,coef,x,param);
end

