function [epipar,obj] = epiapprox_kldiv(pdf_obs,mesh)


mat_inte = [mesh.mlist(2:end)-mesh.mlist(1:end-1),(mesh.mlist(2:end).^2-mesh.mlist(1:end-1).^2)/2];
mat_cont = [[eye(mesh.m-1),zeros(mesh.m-1,1)]+[zeros(mesh.m-1,1),-eye(mesh.m-1)],[diag(mesh.mlist(2:mesh.m)),zeros(mesh.m-1,1)]+[zeros(mesh.m-1,1),-diag(mesh.mlist(2:mesh.m))]];
mat_nonneg = [[eye(mesh.m), diag(mesh.mlist(1:mesh.m))];[zeros(1,mesh.m-1) eye(1) zeros(1,mesh.m-1) mesh.mlist(end)]];

% p=f(x);
% mat_epi=sparse(length(x),2*mesh.m);
% col_index = discretize(x,mesh.mlist);
% row_index = (1:length(x))';
% index=sub2ind(size(mat_epi),row_index,col_index);
% mat_epi(index)=ones(size(x,1),1);
% mat_epi(index+mesh.m*length(x))=x;

n=10;
[quadn,quadw]=lgwt(10,-1,1);
interval=([-eye(mesh.m),zeros(mesh.m,1)]+[zeros(mesh.m,1),eye(mesh.m)])*(mesh.mlist)';
coord = (ones(mesh.m,1)*(flipud(quadn+1)/2)').*((interval)*ones(1,n))+mesh.mlist(1:mesh.m)'*ones(1,n);
weight=repmat(quadw,mesh.m,1)/mesh.m;

mat_epi=zeros(mesh.m*mesh.n,2*mesh.m);
for i=1:mesh.m
    mat_epi((i-1)*mesh.n+1:i*mesh.n,i)=ones(mesh.n,1);
    mat_epi((i-1)*mesh.n+1:i*mesh.n,i+mesh.m)=(coord(i,:))';
end
x=reshape(coord',mesh.m*n,1);
p=pdf_obs(x);

Aeq = [mat_inte;mat_cont];
beq = [1;zeros(mesh.m-1,1)];
A=[-mat_nonneg;mat_nonneg(1,:);mat_nonneg(end,:)];
b=[zeros(mesh.m+1,1);zeros(2,1)];
% options= optimoptions(@lsqlin,'Algorithm','interior-point','MaxIter',100000);
[epipar,obj]=lsqlin(diag(weight)*mat_epi,weight.*p,A,b,Aeq,beq,[],[],[]);
% x0=ones(mesh.m*2,1);
% options= optimoptions(@fmincon,'MaxFunctionEvaluations',100000);
% [epipar,obj]=fmincon(@(par)-weight'*(log(mat_epi*par).*p),x0,A,b,Aeq,beq,[],[],[],options);


epipar=reshape(epipar,mesh.m,length(epipar)/mesh.m);


end