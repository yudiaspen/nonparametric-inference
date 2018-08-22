function epidraw(epipar, mesh, color)
ydim=epipar(:,1)+epipar(:,2).*mesh.mlist(2:end)';
ydim=[epipar(1,1);ydim];
plot(mesh.mlist,ydim,color);
hold on
ylim([0,2.5])
end