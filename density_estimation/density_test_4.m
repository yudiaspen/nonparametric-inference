% Correlated normal
clear

% fast function test
fun_map = @(x) [exp(x(:,1)+x(:,2)+1),(x(:,1)-1).^2,(x(:,2)-1).^2];
num_dim_y_full = 3;num_dim_x_full = 2;

% true density
par_loc = ones(1,num_dim_x_full)*-.1;
par_scale= 0.2*ones(1,num_dim_x_full);
par_corr = 0.5;
par = [par_loc,par_scale,par_corr];
fun_pdf = @(x,par) copula_nataf(x,par(1:num_dim_x_full),par(num_dim_x_full+1:num_dim_x_full*2),par(num_dim_x_full*2+1:end));

% generate data points
rng(1);
obs_xmvn_pool_full = mvnrnd(zeros(1,num_dim_x_full),diag(ones(1,num_dim_x_full).^2)+flip(diag(par_corr*ones(1,num_dim_x_full))),1e6);
obs_x_pool_full = normcdf(obs_xmvn_pool_full);
for i =1:num_dim_x_full, obs_x_pool_full(:,i)=logninv(obs_x_pool_full(:,i),par_loc(i),par_scale(i)); end
obs_y_pool_full = fun_map(obs_x_pool_full);

% naive integration
fun_pdf_smp = @(x) 1/(2^num_dim_x_full);

% generate numerical sampling points
rng(2);
smp_x_full = [rand(1e6,1)*2,rand(1e6,1)*2];
smp_y_full = fun_map(smp_x_full);
if ~exist('pdf_aux','var')
    pdf_aux=ones(size(smp_x_full,1),1);
end

% estimation specification
num_sample=400;
num_obs = 250;
num_smp = 1e5;
num_smp_iter = 1e4;
num_sample_test=2;

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
indices_x=[1,2];num_dim_x=length(indices_x);
density_estimate_p;

function [pdf_joint]=copula_nataf(x,par_loc_x,par_scale_x,par_corr_x)
    norm_x = zeros(size(x));
    for i=1:size(x,2), norm_x(:,i)=cdf2normx(logncdf(x(:,i),par_loc_x(i),par_scale_x(i))); end
    pdf_joint = lognpdf(x(:,1),par_loc_x(1),par_scale_x(1)).*lognpdf(x(:,2),par_loc_x(2),par_scale_x(2))...
                ./ prod(normpdf(norm_x),2) .* mvnpdf(norm_x,zeros(1,size(x,2)),diag(ones(1,size(x,2)))+flip(diag(ones(1,size(x,2))*par_corr_x)));
end