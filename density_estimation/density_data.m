% 1. Data Preparation
mean_y_true = mean(smp_y_full(:,indices_y));

% store normalized data in obs_total_y; obs_pool_y is for convenience
obs_total_y=cell(num_sample,1);obs_pool_y=zeros(num_obs*num_sample,num_dim_y);
for i=1:num_sample
    obs_total_y{i}=obs_y_pool_full((i-1)*num_obs+1:i*num_obs,indices_y)./(ones(num_obs,1)*mean_y_true)-1;
    obs_pool_y((i-1)*num_obs+1:i*num_obs,:)=obs_total_y{i};
end
% smp y mean is set to zero
smp_y=smp_y_full(1:num_smp,indices_y)./(ones(num_smp,1)*mean_y_true(indices_y))-1;