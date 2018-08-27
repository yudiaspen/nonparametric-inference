function [val] = gmmobj_corr(corr,p_mar,x,coef,moment,weight)

num_x=size(x,2);

p=mvnpdf(x,zeros(1,num_x),diag(ones(1,num_x))+flip(diag(ones(1,num_x)*corr))).*p_mar;

diff=(p'*coef)'-moment;

val=diff'*weight*diff;


end
