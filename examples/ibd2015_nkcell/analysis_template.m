% stop the execution temporarily for 1-10 minutes based on the job ID
% this is to avoid license issues when running many jobs in the same time
n_time = abs(mod(case_id, 10))*60;
pause(n_time);

% set up parallel computing environment
pc = parcluster('local');
parpool(pc, str2num(getenv('SLURM_CPUS_ON_NODE')));  %#ok<ST2NM>

% specify file prefix
file_prefix = strcat(gwas_name,'_',net_name,'_',cis_name);

% specify total number of autosome protein-coding genes
data.num_gene = ngene;

% specify GWAS sample size
data.num_sam = nsam;

% turn off RSS program progress bar
options.verbose = false;

% set the convergence tolerance of RSS program
options.tolerance = 1e-3;

% set the maximum wall time for RSS program (unit: seconds)
options.max_walltime = 12*3600;

% create a folder for output files if not exist
out_path = strcat(out_path,file_prefix,'/');
if exist(out_path,'dir') ~= 7
  mkdir(out_path);
end

% add search paths
addpath(strcat(src_path,'rss/src_vb/'));  % RSS-E codes
addpath(strcat(src_path,'rss-net/src/')); % RSS-NET codes

% specify single-SNP GWAS summary statistics, if they are not specified
if ~isfield(data,'sumstats_file')
  data.sumstats_file = strcat(dat_path,gwas_name,'_sumstat.mat');
end

% specify network-related data files, if they are not specified
if ~isfield(data,'snp2net_file')
  data.snp2net_file = strcat(dat_path,gwas_name,'_',net_name,'_snp2net.mat');
end

if ~isfield(data,'snp2gene_file')
  data.snp2gene_file = strcat(dat_path,gwas_name,'_snp2gene.mat');
end

if ~isfield(data,'gene2gene_file') 
  data.gene2gene_file = strcat(dat_path,net_name,'_gene2gene.mat');
end

if ~isfield(data,'snp2gene_cis_file')
  data.snp2gene_cis_file = strcat(dat_path,gwas_name,'_',cis_name,'_snp2gene_cis.mat');
end

% do not use gene expression in current RSS-NET analysis
data.expression_file = [];

% specify network-related parameters, if they are not specified 
if ~isfield(options,'snp2gene_par')
  options.snp2gene_par = 'dist_bin';
end

if ~isfield(options,'snp2gene_opt')
  options.snp2gene_opt = 1e6;
end

if ~isfield(options,'gene2gene_par')
  options.gene2gene_par = 'asis';
end

if ~isfield(options,'gene2gene_opt') 
  options.gene2gene_opt = [];
end

if ~isfield(options,'snp2gene_cis_par')
  options.snp2gene_cis_par = 'asis';
end

if ~isfield(options,'snp2gene_cis_opt')
  options.snp2gene_cis_opt = [];
end

% specify RSS-E whole-genome baseline results
% source: https://doi.org/10.1038/s41467-018-06805-x
base_file = strcat(dat_path,gwas_name,'_null_seed_459_squarem_step2.mat');

% initialize RSS-NET with RSS-E baseline results
rsse_base = matfile(base_file);

logw0_vec  = rsse_base.logw;
alpha0_mat = rsse_base.alpha;
mu0_mat    = rsse_base.mu;

data.num_snp = size(mu0_mat,1); % total # of genome-wide SNPs

[~,mi] = max(logw0_vec(:));
alpha0 = alpha0_mat(:,mi);
mu0    = mu0_mat(:,mi);

if (length(alpha0) ~= data.num_snp) || (length(mu0) ~= data.num_snp)
  error('Inconsistent number of SNPs in baseline results ...');
end

options.alpha = alpha0;
options.mu    = mu0;

clear base_file rsse_base alpha0* mu0* logw0 mi;

% specify hyper-parameters for RSS-NET, using `make_hpgrid.m`
param_data = make_hpgrid(hyper_data);

hyper.theta0 = param_data(case_id,1);
hyper.theta  = param_data(case_id,2);
hyper.eta    = param_data(case_id,3);
hyper.rho    = param_data(case_id,4);

clear param_* hyper_*;

% run RSS-NET analysis
tic;
[logw,alpha,mu,s,param] = rss_net(data,hyper,options);
run_time = toc;

clear data hyper options;

% create name of output file for this analysis
out_file = strcat(out_path,file_prefix,'_out_',num2str(case_id),'.mat');

% save RSS-NET analysis output
theta0 = param(1,1);
theta  = param(1,2);
sigb   = param(1,3);
sige   = param(1,4);

save(out_file,'logw','alpha','mu','s','run_time',...
     'theta0','theta','sigb','sige');

% terminate existing parallel session
delete(gcp('nocreate'));
