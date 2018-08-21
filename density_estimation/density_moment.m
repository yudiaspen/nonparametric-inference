% 2. Generate Moments List
momentlist_full=zeros(1,num_moment);
for i=1:num_moment
momentlist_temp=nmultichoosek(1:num_dim_y,i);
momentlist_full=[momentlist_full;[momentlist_temp,zeros(size(momentlist_temp,1),num_moment-i)]];
end

mat_temp_index=(1:size(momentlist_full,1))'*ones(1,num_moment);
mat_temp=[momentlist_full(:),mat_temp_index(:)];
momentlist=sparse(mat_temp(mat_temp(:,1)~=0,2),mat_temp(mat_temp(:,1)~=0,1),ones(sum(mat_temp(:,1)~=0),1));

% only first two moments
moment_select=filter_moment_select(momentlist);
moment_select(1)=0;

%-----------------------------------------------------------------------------
tic
% generate observation moment
moment_obs_y = zeros(size(momentlist,1),num_sample);
for iteri=1:num_sample
    obs_y=obs_total_y{iteri};
for i=1:size(momentlist,1)
    moment_obs_y(i,iteri)=mean(prod(obs_y.^repmat(momentlist(i,:),[num_obs,1]),2));
end
end
toc
% generate pool observation moment
moment_obs_pool_y = mean(moment_obs_y,2);

