function epidraw(epipar, mesh, color)
epipar=reshape(epipar,mesh.m,length(epipar)/mesh.m);
ydim=epipar(:,1)+epipar(:,2).*mesh.mlist(2:end)';
ydim=[epipar(1,1)+epipar(1,2)*mesh.m0;ydim];
plot(mesh.mlist,ydim,color);
hold on
% ylim([0,2.5])
end