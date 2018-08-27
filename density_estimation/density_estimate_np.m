% moment_obs_y(:,1)=moment_pdfsmp_y;

% mesh_z size
mesh_size = 64;

mesh_x= struct('m',{},'m0',{},'mend',{},'mlist',{});
for i = indices_x_np
    mesh_x(i).m=mesh_size;mesh_x(i).m0=bound_x_min(i);mesh_x(i).mend=bound_x_max(i);
    mesh_x(i).mlist = mesh_x(i).m0:1/mesh_x(i).m*(mesh_x(i).mend-mesh_x(i).m0):mesh_x(i).mend;
end

bound_z_max=sum(max(mat_chol_inv.*[2; 2],mat_chol_inv.*[0; 0]),1);
bound_z_min=sum(min(mat_chol_inv.*[2; 2],mat_chol_inv.*[0; 0]),1);

mesh_z= struct('m',{},'m0',{},'mend',{},'mlist',{});
for i = indices_x_np
    mesh_z(i).m=mesh_size;mesh_z(i).m0=bound_z_min(i);mesh_z(i).mend=bound_z_max(i);
    mesh_z(i).mlist = mesh_z(i).m0:1/mesh_z(i).m*(mesh_z(i).mend-mesh_z(i).m0):mesh_z(i).mend;
end


[quad_z,quad_w]=epi_quad(mesh_z);

grid_dense=combvec(quad_z(:,1)',quad_z(:,2)'); % generate z combinatins
grid_dense_pdf = det(mat_chol)*fun_pdf((grid_dense'*mat_chol),par); % calculate fx using x = Az
quad_p=[reshape(grid_dense_pdf,length(quad_z(:,1)),[])*quad_w(:,2),reshape(grid_dense_pdf,length(quad_z(:,2)),[])'*quad_w(:,1)];

par_epi_z=zeros(mesh_size*2,num_dim_x_full);
for i=indices_x_np,par_epi_z(:,i)=epiapprox(quad_p(:,i) ,quad_z(:,i),quad_w(:,i),mesh_z(i)); end

num_param_np=num_dim_x_np*2*mesh_size; % number of parameters

mesh_z_group = zeros(size(smp_z));
for i = indices_x_np
    mesh_z_group(:,i)=discretize(smp_z(:,i),mesh_z(i).mlist);
end

mesh_x_group = zeros(size(smp_z));
for i = indices_x_np
    mesh_x_group(:,i)=discretize(smp_x_full(1:num_smp,i),mesh_x(i).mlist);
end

fun_pdf_np=@(x,par_epi_z) det(mat_chol_inv)*prod(epipdf(x,par_epi_z,mesh_z_group),2);

pdf_true = zeros(num_smp,num_dim_x_np); 
for i = 1:num_dim_x_np
    pdf_true(:,i) = (par_epi_z(mesh_z_group(:,i)+mesh_size,i).*smp_z(:,i)+par_epi_z(mesh_z_group(:,i),i));
end

[gmmG_true,diffmat]=gmmjacob_epi(pdf_true,coef_cpl(:,moment_select),smp_z(:,indices_x>0),par_epi_z,mesh_z_group);
% [gmmG_true_num,diffmat_num]=gmmjacob_num(fun_pdf_np,coef_cpl(:,moment_select),smp_z(:,indices_x>0),par_epi_z);
cov_true=inv(gmmG_true'*((sigmadata+num_obs/num_smp_iter*sigmasimulation)\gmmG_true))/num_obs;

est_total_np=zeros(num_sample_test,num_param_np);
rng(2)
est_init_np= rand(2*mesh_size,num_dim_x_np)*2-1+par_epi_z;
tic
options = optimoptions('fmincon','TolFun',1e-6,'MaxFunEvals',5e5,'MaxIter',5e4,'GradObj','on','Display','iter');
obj_val = zeros(num_sample_test,1); obj_fit = zeros(num_sample_test,1);%recording, objective function value

mat_inte=cell(num_dim_x_np,1);mat_cont=cell(num_dim_x_np,1);mat_nonneg=cell(num_dim_x_np,1);
for i = indices_x_np
    mat_inte{i} = [mesh_z(i).mlist(2:end)-mesh_z(i).mlist(1:end-1),(mesh_z(i).mlist(2:end).^2-mesh_z(i).mlist(1:end-1).^2)/2];
    mat_cont{i} = [[eye(mesh_size-1),zeros(mesh_size-1,1)]+[zeros(mesh_size-1,1),-eye(mesh_size-1)],[diag(mesh_z(i).mlist(2:mesh_size)),zeros(mesh_size-1,1)]+[zeros(mesh_size-1,1),-diag(mesh_z(i).mlist(2:mesh_size))]];
    mat_nonneg{i} = [[eye(mesh_size), diag(mesh_z(i).mlist(1:mesh_size))];[zeros(1,mesh_size-1) eye(1) zeros(1,mesh_size-1) mesh_z(i).mlist(end)]];
end
Aeq = [blkdiag(mat_inte{:});blkdiag(mat_cont{:})];
beq = [ones(num_dim_x_np,1);zeros((mesh_size-1)*num_dim_x_np,1)];
A=-blkdiag(mat_nonneg{:});
b=zeros((mesh_size+1)*num_dim_x_np,1);



for j=1:num_sample_test
    rng(j)
    rand_sim= randsample(num_smp,num_smp_iter);
    gmmweight=inv(sigmadata+num_obs/length(rand_sim)*sigmasimulation);
    fitness_np=@(param)gmmobj_epi(param,smp_z(rand_sim,indices_x),coef_cpl(rand_sim,moment_select)*num_smp/num_smp_iter,moment_obs_y(moment_select,j),gmmweight,mesh_z_group(rand_sim,:));
    fitness=@(param)gmmobj(param,smp_z(rand_sim,indices_x),coef_cpl(rand_sim,moment_select)*num_smp/num_smp_iter,moment_obs_y(moment_select,j),gmmweight,@(x,par_epi_z) prod(epipdf(x,par_epi_z,mesh_z_group(rand_sim,:))),2);
    
%     fitness_np=fitness;
    
    tic
    [est_z]=  fmincon(fitness_np,est_init_np,A,b,Aeq,beq ,par_epi_z-1,par_epi_z+1,[],options);
    est_total_np(i,:)=reshape(est_z,num_param_np,1);
    obj_val(i) = fitness_np(est_z);obj_fit(i) = fitness_np(par_epi_z);
    toc
   
    
    [quad_x,quad_w]=epi_quad(mesh_x);
    
    quad_p=epiconv(quad_x,par_epi_z,mesh_z,mat_chol);
    
    grid_dense=combvec(quad_x(:,1)',quad_x(:,2)'); % generate z combinatins
    mesh_z_group_temp = zeros(size(grid_dense'));temp_z=grid_dense'/mat_chol;
    for i = indices_x_np
        mesh_z_group_temp(:,i)=discretize(temp_z(:,i),mesh_z(i).mlist);
    end
    grid_dense_pdf = det(mat_chol_inv)*prod(epipdf(temp_z,par_epi_z,mesh_z_group_temp),2); % calculate fx using x = Az
    quad_p=[reshape(grid_dense_pdf,length(quad_x(:,1)),[])*quad_w(:,2),reshape(grid_dense_pdf,length(quad_x(:,2)),[])'*quad_w(:,1)];
    
    
    est_x = zeros(mesh_size*2,length(indices_x_np));
    for i=indices_x_np,est_x(:,i)=epiapprox(quad_p(:,i) ,quad_x(:,i),quad_w(:,i),mesh_x(i)); end
    
   
    [p_mar]=epipdf(smp_x_full(rand_sim,:),est_x,mesh_x_group(rand_sim,:));
    xcdf = epicdf(smp_x_full(rand_sim,:),est_x,mesh_x,mesh_x_group(rand_sim,:));
    norm_x=cdf2normx(xcdf);
    
%     p_mar = [lognpdf(smp_x_full(rand_sim,1),par_loc(1),par_scale(1)),lognpdf(smp_x_full(rand_sim,2),par_loc(2),par_scale(2))];
%     xcdf =[logncdf(smp_x_full(rand_sim,1),par_loc(1),par_scale(1)),logncdf(smp_x_full(rand_sim,2),par_loc(2),par_scale(2))];
%     norm_x=cdf2normx(xcdf);
    
    p_mar=prod(p_mar,2)./ prod(normpdf(norm_x),2);
    
    est_init_corr=0;
    fitness_corr=@(corr)gmmobj_corr(corr,p_mar,norm_x,coef_cpl(rand_sim,moment_select)*num_smp/num_smp_iter,moment_obs_y(moment_select,j),gmmweight);
    [est_corr]=  fmincon(fitness_corr,est_init_corr,[],[],[],[] ,-.9,.9,[]);
    
end

% plot z: b true g approx r estimate
epidraw(par_epi_z(:,1),mesh_z(1),'g');epidraw(par_epi_z(:,2),mesh_z(2),'g');
% epidraw(est_total_np(i,mesh_size*2*0+1:mesh_size*1*2),mesh_z(1),'r');epidraw(est_total_np(i,mesh_size*2*1+1:mesh_size*2*2),mesh_z(2),'r');
plot(-1:0.01:2.5,pdf(fitdist(obs_z_pool_full(1:1e6,1),'kernel'),-1:0.01:2.5),'b');hold on;plot(-1:0.01:2.5,pdf(fitdist(obs_z_pool_full(1:1e6,2),'kernel'),-1:0.01:2.5),'b');

% plot x: b true g approx r estimate
epidraw(est_x(:,1),mesh_x(1),'r');epidraw(est_x(:,2),mesh_x(2),'r');
plot((0:.01:2)',fun_mar{1}((0:.01:2)),'b');hold on;plot((0:.01:2)',fun_mar{2}((0:.01:2)),'b');
p=epiconv([(0:.01:2)',(0:.01:2)'],par_epi_z,mesh_z,mat_chol);hold on;plot((0:.01:2)',p(:,1),'c');hold on;plot((0:.01:2)',p(:,2),'c');