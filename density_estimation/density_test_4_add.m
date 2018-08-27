% Correlated normal
clear

% fast function test
fun_map = @(x) [exp(x(:,1)+x(:,2)+1),(x(:,1)-1).^2,(x(:,2)-1).^2];
num_dim_y_full = 3;num_dim_x_full = 2;

% true density
par_loc = [-.25,-.25];
par_scale= [.25,.25];
par_corr = 0;
par = [par_loc,par_scale,par_corr];
mat_var = diag(ones(1,num_dim_x_full).^2)+flip(diag(par_corr*ones(1,num_dim_x_full)));
% fun_pdf = @(x,par) copula_nataf(x,par(1:num_dim_x_full),par(num_dim_x_full+1:num_dim_x_full*2),par(num_dim_x_full*2+1:end));
fun_pdf = @(x,par) copulapdf('Gaussian',[logncdf(x(:,1),par(1),par(3)),logncdf(x(:,2),par(2),par(4))],par(5)).*lognpdf(x(:,1),par(1),par(3)).*lognpdf(x(:,2),par(2),par(4));
fun_mar = {@(x)lognpdf(x,par_loc(1),par_scale(1)),@(x)lognpdf(x,par_loc(2),par_scale(2))};

% generate data points
rng(1);
obs_xmvn_pool_full = mvnrnd(zeros(1,num_dim_x_full),diag(ones(1,num_dim_x_full).^2)+flip(diag(par_corr*ones(1,num_dim_x_full))),1e6);
obs_x_pool_full = normcdf(obs_xmvn_pool_full);
for i =1:num_dim_x_full, obs_x_pool_full(:,i)=logninv(obs_x_pool_full(:,i),par_loc(i),par_scale(i)); end
obs_y_pool_full = fun_map(obs_x_pool_full);
 %mat_chol = chol(diag(ones(1,num_dim_x_full).^2)+flip(diag(par_corr*ones(1,num_dim_x_full))));
mat_chol = chol(corr(obs_x_pool_full));
obs_z_pool_full=obs_x_pool_full/mat_chol;

% obs_copula=copularnd('Gaussian',par_corr,1e6);
% for i =1:num_dim_x_full, obs_x_pool_full(:,i)=logninv(obs_copula(:,i),par_loc(i),par_scale(i));end
% obs_z_pool_full=obs_x_pool_full/mat_chol;  
% pd1=fitdist(obs_z_pool_full(1:1e5,1),'kernel');
% pd2=fitdist(obs_z_pool_full(1:1e5,2),'kernel');
% pdf1=fun_pdf(obs_x_pool_full(1:1e5,:),par);
% pdf2=pdf(pd1,obs_z_pool_full(1:1e5,1)).*pdf(pd2,obs_z_pool_full(1:1e5,2))/det(mat_chol);
% scatter(pdf1,pdf2,.01);

% naive integration
bound_x_max = ones(1,num_dim_x_full)*2;
bound_x_min = ones(1,num_dim_x_full)*0;
fun_pdf_smp = @(x) 1/(prod(bound_x_max-bound_x_min));

% generate numerical sampling points
rng(2);
smp_x_full = rand(1e6,num_dim_x_full).*(bound_x_max-bound_x_min)+ones(1e6,1)*bound_x_min;
smp_y_full = fun_map(smp_x_full);
if ~exist('pdf_aux','var')
    pdf_aux=ones(size(smp_x_full,1),1);
end

% estimation specification
num_sample=200;
num_obs = 2500;
num_smp = 1e5;
num_smp_iter = 1e4;
num_sample_test=1;

% only consider observed y
indices_y=[1,2,3];num_dim_y=length(indices_y); 
density_data;

% only consider selected moments
num_moment=4;
filter_moment_select= @(my)(sum(my,2)<=10);
density_moment;

% coefficient computation
density_coef;

% only consider estimated x
% moment_obs_y(:,1)=moment_pdfsmp_y;
indices_x=[1,2];num_dim_x=length(indices_x);
density_estimate_p;

% nonparametric estimation
indices_x_np=[1,2];


mat_chol_inv=inv(mat_chol);
smp_z = smp_x_full(1:num_smp,:)/mat_chol;

num_dim_x_np=numel(indices_x_np); % x dimension



density_estimate_np;

function [pdf_joint]=copula_nataf(x,par_loc_x,par_scale_x,par_corr_x)
    norm_x = zeros(size(x));
    for i=1:size(x,2), norm_x(:,i)=cdf2normx(normcdf(x(:,i),par_loc_x(i),par_scale_x(i))); end
    pdf_joint = lognpdf(x(:,1),par_loc_x(1),par_scale_x(1)).*lognpdf(x(:,2),par_loc_x(2),par_scale_x(2))...
                ./ prod(normpdf(norm_x),2) .* mvnpdf(norm_x,zeros(1,size(x,2)),diag(ones(1,size(x,2)))+flip(diag(ones(1,size(x,2))*par_corr_x)));
end