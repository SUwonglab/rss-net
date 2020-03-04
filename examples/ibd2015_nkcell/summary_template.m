% add search path
addpath(strcat(src_path,'rss/src_vb/')); % use compute_pip.m
addpath(strcat(src_path,'rss-net/src/')); % use calc_bf.m

% specify file prefix: [gwas]_[net]_[cis]
out_prefix = strcat(gwas_name,'_',net_name,'_',cis_name);

% aggregate results across hyper-parameters
num_par = length(par_set);

theta0 = zeros(num_par,1);
theta  = zeros(num_par,1);
sigb   = zeros(num_par,1);
sige   = zeros(num_par,1);
logw   = zeros(num_par,1);
time   = zeros(num_par,1);

for j=1:num_par
  % load full results
  res = matfile(strcat(res_path,out_prefix,'_out_',num2str(par_set(j)),'.mat'));

  % extract model-level results
  theta0(j) = res.theta0;
  theta(j)  = res.theta;
  sigb(j)   = res.sigb;
  sige(j)   = res.sige;
  logw(j)   = res.logw;
  time(j)   = res.run_time;

  % preallocate variational output
  if j==1
    logw_net  = zeros(num_par,1);
    alpha_net = zeros(length(res.alpha),num_par);
  end

  % extract variational estimates
  logw_net(j)    = res.logw;
  alpha_net(:,j) = res.alpha;

  fprintf('Results of file %s are extracted ...\n',num2str(par_set(j)));
end

% compute enrichment BFs
param = [theta0 theta sigb sige];

[log10_bf,log10_bf_nt,log10_bf_ns,log10_bf_ts] = calc_bf(logw,param);
clear param;

% specify genomic locus info
gene_seg_path = strcat(dat_path,gwas_name,'_gene_grch37.mat');

% summarize locus-level results (gene +/- 100kb)
P1_gene = compute_pip(gene_seg_path, logw_net, alpha_net);
clear logw_net alpha_net;

gene_seginfo = matfile(gene_seg_path);

gene_chr   = gene_seginfo.segchr;
gene_start = gene_seginfo.segstart;
gene_stop  = gene_seginfo.segstop;

% save results
out_file = strcat(out_prefix,'_results_model.mat');
save(out_file,'log10_bf','log10_bf_nt','log10_bf_ns','log10_bf_ts',...
     'theta0','theta','sigb','sige','logw','time');

out_file = strcat(out_prefix,'_results_gene.mat');
save(out_file,'gene_chr','gene_start','gene_stop','P1_gene');

