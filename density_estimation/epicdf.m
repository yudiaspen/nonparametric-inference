function [cd]=epicdf(x,param,mesh,g)

num_x=size(x,2);
msize=size(param,1)/2;

parb=param(1:msize,:);
para=param(msize+1:end,:);

mlist=reshape([mesh(:).mlist],[],num_x);
mright=mlist(2:end,:);
mleft=mlist(1:end-1,:);

mleft2=mleft.^2;

cd_cum=[zeros(1,num_x);cumsum(para.*(mright.^2-mleft2)/2+parb.*(mright-mleft));];

cd=zeros(size(x));

for i=1:num_x
cd(:,i) = cd_cum(g(:,i),i)+ para(g(:,i),i).*(x(:,i).^2-mleft2(g(:,i),i))/2+parb(g(:,i),i).*(x(:,i)-mleft(g(:,i),i));
end

end

