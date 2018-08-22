% 3. Coefficient Computation
% integrate gx using monte carlo integration
weight=1./fun_pdf_smp(smp_x_full)/num_smp./pdf_aux(1:num_smp,:);
coefpure=zeros(num_smp,size(momentlist,1));
for i=1:size(momentlist,1)

    coefpure(:,i)=prod(smp_y(1:num_smp,:).^repmat(momentlist(i,:),[num_smp,1]),2).*weight;

end 

coefpure=coefpure.*(pdf_aux(1:num_smp,:)*ones(1,size(momentlist,1)));
%-----------------------------------------------------------------------------
tic
% generate integrated population moment
moment_pdfsmp_y = zeros(size(momentlist,1),1);
pdf_smp = fun_pdf(smp_x_full(1:num_smp,:),par);


% sample y pdf
for i=1:size(momentlist,1)
    moment_pdfsmp_y(i)=dot(prod(smp_y(1:num_smp,:).^repmat(momentlist(i,:),[num_smp,1]),2).*weight,pdf_smp.*pdf_aux(1:num_smp,:));
end
toc

%------------------------------------------------------------------------------
tic
 [sigmadata]=weight_obs(obs_pool_y(1:1e4,:),moment_obs_pool_y(moment_select),momentlist(moment_select,:));
toc

sigmasimulation=weight_smp(smp_y(1:1e4,:),1./fun_pdf_smp(smp_x_full(1:1e4,:)).*fun_pdf(smp_x_full(1:1e4,:),par),moment_obs_pool_y(moment_select),momentlist(moment_select,:));