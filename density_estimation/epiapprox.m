function [epipar,obj] = epiapprox(p,x,w,mesh)


mat_inte = [mesh.mlist(2:end)-mesh.mlist(1:end-1),(mesh.mlist(2:end).^2-mesh.mlist(1:end-1).^2)/2];
mat_cont = [[eye(mesh.m-1),zeros(mesh.m-1,1)]+[zeros(mesh.m-1,1),-eye(mesh.m-1)],[diag(mesh.mlist(2:mesh.m)),zeros(mesh.m-1,1)]+[zeros(mesh.m-1,1),-diag(mesh.mlist(2:mesh.m))]];
mat_nonneg = [[eye(mesh.m), diag(mesh.mlist(1:mesh.m))];[zeros(1,mesh.m-1) eye(1) zeros(1,mesh.m-1) mesh.mlist(end)]];
Aeq = [mat_inte;mat_cont];
beq = [1;zeros(mesh.m-1,1)];
A=[-mat_nonneg;mat_nonneg(1,:);mat_nonneg(end,:)];
b=[zeros(mesh.m+1,1);zeros(2,1)];

% p=f(x);
% mat_epi=sparse(length(x),2*mesh.m);
% col_index = discretize(x,mesh.mlist);
% row_index = (1:length(x))';
% index=sub2ind(size(mat_epi),row_index,col_index);
% mat_epi(index)=ones(size(x,1),1);
% mat_epi(index+mesh.m*length(x))=x;

n=size(x,1)/mesh.m;
mat_epi=zeros(mesh.m*n,2*mesh.m);
for i=1:mesh.m
    mat_epi((i-1)*n+1:i*n,i)=ones(n,1);
    mat_epi((i-1)*n+1:i*n,i+mesh.m)=x((i-1)*n+1:i*n)';
end




% options= optimoptions(@lsqlin,'Algorithm','interior-point','MaxIter',100000);
[epipar,obj]=lsqlin(diag(w)*mat_epi,w.*p,A,b,Aeq,beq,[],[],[]);
% x0=ones(mesh.m*2,1);
% options= optimoptions(@fmincon,'MaxFunctionEvaluations',100000);
% [epipar,obj]=fmincon(@(par)-weight'*(log(mat_epi*par).*p),x0,A,b,Aeq,beq,[],[],[],options);


% epipar=reshape(epipar,mesh.m,length(epipar)/mesh.m);


end