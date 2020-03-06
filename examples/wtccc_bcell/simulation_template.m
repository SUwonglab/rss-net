% add search paths
addpath(strcat(src_path,'rss/src_vb/'));  % RSS-E codes
addpath(strcat(src_path,'rss-net/src/')); % RSS-NET codes

% specify GWAS and network info
gwas_name = 'WTCCC';
net_name  = 'Primary_B_cells_from_peripheral_blood';

% get the number of individuals and SNPs
data.num_sam = 1458;
data.num_snp = 348965;

% specify the number of genes
data.num_gene = 18334;

% specify single-SNP GWAS summary statistics
data.sumstats_file = strcat(dat_path,trial_name,'_data.mat');

% specify input network data files
data.snp2net_file      = strcat(dat_path,gwas_name,'_',net_name,'_snp2net.mat');
data.snp2gene_file     = strcat(dat_path,gwas_name,'_snp2gene.mat');
data.gene2gene_file    = strcat(dat_path,net_name,'_gene2gene.mat');
data.snp2gene_cis_file = strcat(dat_path,gwas_name,'_',net_name,'_snp2gene_cis.mat');
data.expression_file   = [];

% specify genomic locus information
seg_path = strcat(dat_path,'wtccc_tidy_snp2gene_ncbi35togrch37.mat');

% specify hyper-parameter grid for RSS-NET analysis
hyper.theta0 = true_theta0;
hyper.sigma0 = true_sigma0;

if (true_theta==0)
  hyper.theta = (0:0.25:1);
else
  hyper.theta = [0 ((true_theta-0.5):0.25:(true_theta+0.5))];
end

if (true_sigma==0)
  hyper.sigma = (0:0.25:1);
else
  hyper.sigma = [0 ((true_sigma-0.5):0.25:(true_sigma+0.5))];
end

% specify parameters for creating network matrix 
options.snp2gene_par = 'dist_bin';
options.snp2gene_opt = 1e6;

options.gene2gene_par = 'asis';
options.gene2gene_opt = [];

options.snp2gene_cis_par = 'asis';
options.snp2gene_cis_opt = [];

% set initial variational parameters to ensure reproducibility
rng(seed, 'twister');
alpha0 = rand(data.num_snp,1);
alpha0 = alpha0 ./ sum(alpha0);

rng(seed+1, 'twister');
mu0 = randn(data.num_snp,1);

options.alpha = alpha0;
options.mu    = mu0;

clear alpha0 mu0 seed;

% turn off the display of RSS-NET VB progress
options.verbose = false;

% run RSS-NET analysis
tic;
[logw,alpha,~,~,param] = rss_net(data,hyper,options);
run_time = toc;

clear data hyper options;

% assess network-level enrichments
log10_bf = calc_bf(logw,param);

% assess gene-level associations
pp = compute_pip(seg_path,logw,alpha);

clear logw alpha mu s param;

% specify output file name
rss_out = strcat(out_path,trial_name,'_results.mat');

% save output
save(rss_out,'run_time','log10_bf','pp'); 
