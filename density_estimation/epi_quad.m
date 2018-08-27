function [x,w] = epi_quad(mesh)
n=10;
dim=length(mesh);
w=zeros(mesh(1).m*n,dim);x=zeros(mesh(1).m*n,dim);
for i=1:length(mesh)

[quadn,quadw]=lgwt(n,-1,1);
interval=([-eye(mesh(i).m),zeros(mesh(i).m,1)]+[zeros(mesh(i).m,1),eye(mesh(i).m)])*(mesh(i).mlist)';
coord = (ones(mesh(i).m,1)*(flipud(quadn+1)/2)').*((interval)*ones(1,n))+mesh(i).mlist(1:mesh(i).m)'*ones(1,n);
w(:,i)=reshape(quadw*(mesh(i).mlist(2:end)-mesh(i).mlist(1:end-1))/2,[],1);
x(:,i)=reshape(coord',mesh(i).m*n,1);
end

end