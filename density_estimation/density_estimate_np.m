
num_param_np=num_dim_x*2*mesh.m+num_corr; % number of parameters


[gmmG_true,diffmat]=gmmjacob(coef_cpl(:,moment_select),smp_x_full(1:num_smp,indices_x>0),par_epi);
cov_true=inv(gmmG_true'*((sigmadata+0*num_obs/num_smp_iter*sigmasimulation)\gmmG_true))/num_obs;

est_total=zeros(num_sample_test,num_param);

tic
options = optimoptions('fmincon','TolFun',1e-12,'MaxFunEvals',5e4,'MaxIter',1e4,'GradObj','off','Display','iter');
est_init=[ones(1,num_dim_x)*0.1,ones(1,num_dim_x)*.1,ones(1,num_corr)*0];
obj_val = zeros(num_sample_test,1); %recording, objective function value
for i=1:num_sample_test
    rng(i)
    rand_sim= randsample(num_smp,num_smp_iter);
    gmmweight=inv(sigmadata+num_obs/length(rand_sim)*sigmasimulation);
    fitness=@(param)gmmobj(param,smp_x_full(rand_sim,indices_x>0),coef_cpl(rand_sim,moment_select)*num_smp/num_smp_iter,moment_obs_y(moment_select,i),gmmweight,fun_pdf);
    tic
    [est_total_temp]=  fmincon(fitness,est_init,[],[],[],[] ,[par_loc_x-.5,par_scale_x*.1,ones(1,num_corr)*-.9],[par_loc_x+.5,par_scale_x*10,ones(1,num_corr)*.9],[],options);
    est_total(i,:)=est_total_temp;
    obj_val(i) = fitness(est_total_temp);
    toc
end
toc
