function [val,grad] = gmmobj(param,x,coef,moment,weight,fun_pdf)

num_dim_x=size(x,2);
diff=(fun_pdf(x,param(1:num_dim_x),param(num_dim_x+1:2*num_dim_x),param(2*num_dim_x+1:end))'*coef)'-moment;
% Braess
% diff=((ones(size(x))*1/param.*(x<=param))'*coef)'-moment;
val=diff'*weight*diff;

grad=2*diff'*weight*gmmjacob(coef,x,param);
end

