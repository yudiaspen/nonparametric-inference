function p=epiconv(x,param,mesh,c)
c=c';
num_x=size(x,2);
msize = length(mesh(1).mlist)-1;
num_smp=size(x,1);
p=zeros(size(x));
for i=1:num_x
    num_ind=numel(nonzeros(c(i,:)));
    if num_ind>1
        
        
        mlist=reshape([mesh(:).mlist],[],num_ind);
        mlist=mlist.*c(i,:);
        mright=mlist(2:end,:);
        mleft=mlist(1:end-1,:);
        
        
        
        xmat = x(:,i)*ones(1,msize^num_ind);
        [rightvec1,rightvec2]=meshgrid(mright(:,1),mright(:,2));
        rightmat1 = ones(num_smp,1)*reshape(rightvec1,[],1)';
        rightmat2 = ones(num_smp,1)*reshape(rightvec2,[],1)';
        [leftvec1,leftvec2]=meshgrid(mleft(:,1),mleft(:,2));
        leftmat1 = ones(num_smp,1)*reshape(leftvec1,[],1)';
        leftmat2 = ones(num_smp,1)*reshape(leftvec2,[],1)';
        ub=min(rightmat1,xmat-leftmat2);
        lb=max(leftmat1,xmat-rightmat2);
        one_slot = sparse(ub>lb);
        ub(~one_slot)=0;
        lb(~one_slot)=0;
        
        
        parb=param(1:msize,:)./(c(i,:));
        para=param(msize+1:end,:)./(c(i,:).^2);
        [paramat1,paramat2]=meshgrid(para(:,1),para(:,2));
        paravec1 = reshape(paramat1,[],1)';
        paravec2 = reshape(paramat2,[],1)';
        [parbmat1,parbmat2]=meshgrid(parb(:,1),parb(:,2));
        parbvec1 = reshape(parbmat1,[],1)';
        parbvec2 = reshape(parbmat2,[],1)';
        
        term3=(ub.^3-lb.^3)/3.*(-paravec1.*paravec2);
        term2=(ub.^2-lb.^2)/2.*(paravec1.*parbvec2-parbvec1.*paravec2)+(ub.^2-lb.^2)/2.*xmat.*paravec1.*paravec2;
        term1=(ub-lb).*(parbvec1.*parbvec2)+(ub-lb).*xmat.*(paravec2.*parbvec1);
        
        p(:,i)=sum(term3+term2+term1,2);
        
    else
        ind=find(c(i,:)~=0);
        g =discretize(x(:,i),mesh(ind).mlist);
        p(:,i) = param(g+msize,ind)./(c(i,ind).^2).*x(:,i)+param(g,ind)./(c(i,ind));
    end
end

end