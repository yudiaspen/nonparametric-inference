num_dim_x=length(indices_x);
num_corr = length(par_corr);
num_param=num_dim_x*2+num_corr; % number of parameters
if isempty(setdiff(indices_x,1:num_dim_x_full)), coef_cpl=coefpure; else, error('provide conditional density'); end

if isempty(setdiff(indices_x,1:num_dim_x_full)), coef_cpl=coefpure;
else
    indices_x_not=setdiff(indices_x,1:num_dim_x_full);
    coef_cpl = coefpure.*(mvnpdf(smp_x_full(1:num_smp,indices_x_not),par_loc_x(indices_x_not),par_scale_x(indices_x_not))*ones(1,length(moment_select)));
end

% [gmmG_true,diffmat]=gmmjacob(coef_cpl(:,moment_select),smp_x_full(1:num_smp,indices_x),[par_loc(indices_x),par_scale(indices_x)]);
[gmmG_true,diffmat]=gmmjacob_num(fun_pdf,coef_cpl(:,moment_select),smp_x_full(1:num_smp,indices_x),[par_loc(indices_x),par_scale(indices_x),par_corr]);

cov_true=inv(gmmG_true'*((sigmadata+num_obs/num_smp_iter*sigmasimulation)\gmmG_true))/num_obs;

est_total=zeros(num_sample_test,num_param);


options = optimoptions('fmincon','TolFun',1e-12,'MaxFunEvals',5e4,'MaxIter',1e4,'GradObj','off','Display','notify');
% est_init=[ones(1,num_dim_x).*(par(1:num_dim_x)-.25),ones(1,num_dim_x).*(par(num_dim_x+1:2*num_dim_x)*3),ones(1,num_corr)*0];
obj_val = zeros(num_sample_test,1); obj_fit = zeros(num_sample_test,1);%recording, objective function value
for i=1:num_sample_test
    rng(i)
    rand_sim= randsample(num_smp,num_smp_iter);
    gmmweight=inv(sigmadata+num_obs/length(rand_sim)*sigmasimulation);
    fitness=@(param)gmmobj(param,smp_x_full(rand_sim,indices_x),coef_cpl(rand_sim,moment_select)*num_smp/num_smp_iter,moment_obs_y(moment_select,i),gmmweight,fun_pdf);
    
    disp("Sample "+num2str(i));
    flag = 1;
    while flag
        disp("Trial "+num2str(flag));
        est_init = rand(size(par)).*[ones(1,num_dim_x),par_scale*4.8,ones(1,num_corr)*1.9]+[par_loc-.5,par_scale*.1,ones(1,num_corr)*-.95];
        tic
        [est_total_temp]=  fmincon(fitness,est_init,[],[],[],[] ,[par_loc-.5,par_scale*.2,ones(1,num_corr)*-.9],[par_loc+.5,par_scale*10,ones(1,num_corr)*.9],[],options);
        toc
        if fitness(est_total_temp)<fitness(par)+1e-6 || flag>0, flag = 0; else, flag=flag+1; end
    end
    est_total(i,:)=est_total_temp;
    obj_val(i) = fitness(est_total_temp);obj_fit(i) = fitness(par);
    toc
end

disp(strcat('successfully solved:'," ",num2str(sum(obj_val<obj_fit))," ",'of'," ",num2str(num_sample_test)));
cov_est=cov(est_total(obj_val<obj_fit,:));



