function [G,diff_mat]=gmmjacob_epi(p,coef,x,param,g)

diff_mat=zeros(size(x,1),numel(param));

msize=size(param,1)/2;

num_smp=size(x,1);

p_joint=prod(p,2);
for i=1:size(x,2)

    mat_b = zeros(num_smp,msize);
    ind_b = sub2ind([num_smp,msize],1:num_smp,g(:,i)');
    mat_b(ind_b) = ones(num_smp,1);
    
    mat_a = zeros(num_smp,msize);
    ind_a = sub2ind([num_smp,msize],1:num_smp,g(:,i)');
    mat_a(ind_a) = x(:,i);    
    
    p_cond=p_joint./p(:,i);
    
    diff_mat(:,(i-1)*msize*2+1:i*msize*2)=p_cond.*[mat_b,mat_a];
end

G=coef'*diff_mat;

end