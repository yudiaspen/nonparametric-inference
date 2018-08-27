function [val,grad] = gmmobj_epi(param,x,coef,moment,weight,g)

% g=g(1:size(x,1),:);

msize=size(param,1)/2;

% p = zeros(size(x)); 
% 
% for i = 1:size(x,2)
%     p(:,i) = (param(g(:,i)+msize,i).*x(:,i)+param(g(:,i),i));
% end

[p]=epipdf(x,param,g);

diff=(prod(p,2)'*coef)'-moment;

val=diff'*weight*diff;

grad=2*diff'*weight*gmmjacob_epi(p,coef,x,param,g);

end
