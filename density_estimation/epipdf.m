function [p]=epipdf(x,param,g)

msize=size(param,1)/2;
g=g(1:size(x,1),:);
p = zeros(size(x)); 

for i = 1:size(x,2)
    p(:,i) = (param(g(:,i)+msize,i).*x(:,i)+param(g(:,i),i));
end
% p(p<0)=0;

end

