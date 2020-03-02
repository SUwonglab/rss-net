function [logw,alpha,mu,s,param] = rss_net(data,hyper,options)
% USAGE: fit RSS-NET model for a given set of hyperparameters
% INPUT:
%	data: paths of input data files for RSS-NET, structure
%	hyper: hyper-parameter values for RSS-NET, structure
%	options: user-specified behaviour of RSS-NET, structure
% OUTPUT:
%	logw: variational lower bound (up to a constant), num_hyper by 1
%	alpha: variational posterior inclusion probabilities, num_snp_fitted by num_hyper
%	mu: variational posterior means of additive effects, num_snp_fitted by num_hyper
%	s: variational posterior variances of additive effects, num_snp_fitted by num_hyper
%	param: grid of hyper-parameters [theta0 theta sigma0 sigma], num_hyper by 4 

  % If the input variable `options` is not provided, then use default
  % in make_{netmat,expr,logodds}.m and rss_varbvsr_bigmem_squarem.m.
  if ~exist('options','var')
    options = [];
    disp('Default options for RSS-NET are used here ...');
  end

  % Load input data.
  if ~exist('data','var')
    error('Input data `data` must be provided ...');
  end

  % Load hyper-parameters.
  if ~exist('hyper','var')
    error('Hyper-parameters `hyper` must be provided ...');
  end

  % Specify the file of GWAS summary statistics and LD estimates.
  sumstats = data.sumstats_file;

  % Determine whether whole-genome data are fed to RSS-NET.
  % If yes, RSS-NET output can be directly used for further inferences.
  % If no, RSS-NET output should be adjusted before further inferences.

  % For more details on this post-analysis adjustment, please see:
  % 1. https://doi.org/10.1038/s41467-018-06805-x; 
  % 2. https://doi.org/10.1371/journal.pgen.1003770.

  % If the GWAS summary data mat file contains a variable named `snps`,
  % then RSS-NET assumes this mat file does not contain the whole-genome data.
  % To avoid error please do not include `snps` in any whole-genome data.

  % The variable `data.num_snp_fitted` denotes the number of SNPs fitted in RSS-NET.
  % The variable `data.num_snp` denotes the total number of whole-genome SNPs.
  % These two variables are the same if and only if `whole_genome == true`.

  % Finally I applied RSS-NET only to whole-genome summary data to generate
  % simulation and data analysis results for RSS-NET manuscript, thanks to the
  % high-performance computing facilities at University of Chicago (https://rcc.uchicago.edu/)
  % and Stanford University (https://srcc.stanford.edu/).

  sumstats_info = who('-file',sumstats);
  sumstats_data = matfile(sumstats);

  if ~isfield(data,'num_snp_fitted')
    data.num_snp_fitted = data.num_snp;
  end

  % Scenario 1: RSS-NET is applied to genome-wide SNPs.
  whole_genome = ~ismember('snps',sumstats_info);
  whole_genome = whole_genome && (data.num_snp_fitted == data.num_snp);

  if whole_genome
    fprintf('This RSS-NET analysis uses %d genome-wide SNPs ...\n',data.num_snp_fitted);
  end

  % Scenario 2: RSS-NET is applied to near-network SNPs. 
  subset_genome = ismember('snps',sumstats_info);
  subset_genome = subset_genome && (data.num_snp_fitted ~= data.num_snp);

  if subset_genome
    snps = sumstats_data.snps;

    if (length(sumstats_data.snps) ~= data.num_snp_fitted)
      error('Inconsistent number of near-network SNPs ...\n');
    end

    fprintf('This RSS-NET analysis uses %d near-network SNPs ...\n',data.num_snp_fitted);
  end

  % Stop RSS-NET if neither Scenario 1 nor 2 is correctly specified.
  if (~whole_genome) && (~subset_genome)
    error('Invalid SNP specification of RSS-NET ...\n');
  end

  clear sumstats_data sumstats_info whole_genome;

  % Create whole-genome SNP-gene network information matrix.
  net_mat = make_netmat(data,options);

  % Subset network matrix if whole-genome data are not used.
  if subset_genome
    net_mat = net_mat(snps, :);
    clear snps;
  end

  % Obtain normalized context-specific expression values.
  expr_nval = make_expr(data,options); 

  % Compute SNP-gene scaling factors.
  gene_var = full((net_mat.^2) * expr_nval);
  clear net_mat expr_nval;

  % Create a 4-column matrix for all sets of hyper-parameters.
  % a) param == [theta0 theta eta rho] if induced_flag == true
  % b) param == [theta0 theta sigma0 sigma] if induced_flag == false

  [param,induced_flag] = make_hpgrid(hyper);
  clear hyper;
  
  % Preallocate variational output.
  num_hyper      = size(param,1);
  num_snp_fitted = data.num_snp_fitted;
 
  logw  = zeros(num_hyper,1);
  alpha = zeros(num_snp_fitted,num_hyper);
  mu    = zeros(num_snp_fitted,num_hyper);
  s     = zeros(num_snp_fitted,num_hyper);

  % Run variational approximation for each set of hyper-parameters.
  for k=1:num_hyper

    % Specify hyper-parameters related to log odds.
    theta0  = param(k,1);
    theta   = param(k,2);
    logodds = make_logodds(data.snp2net_file,theta0,theta,options);

    % Specify hyper-parameters related to effect sizes.

    % Note that one can either directly specify prior for (sigma0, sigma),
    % or induce this prior from the prior specified for (eta, rho).
    % However, inducing prior for (sigma0, sigma) from (eta, rho) is valid
    % only if the whole-genome data are fed to RSS-NET. 

    if induced_flag && subset_genome
      error('Induced prior can only be applied to whole-genome data ...');
    end

    if induced_flag
      eta = param(k,3);
      rho = param(k,4);

      [sigma0,sigma] = make_sigma(data,gene_var,logodds,eta,rho);

      param(k,3) = sigma0;
      param(k,4) = sigma;
    else
      sigma0 = param(k,3);
      sigma  = param(k,4);
    end
   
    % Compute prior SD for each SNP. 
    sigbeta = sqrt(sigma0^2 + sigma^2*gene_var);  

    % Display the current status.
    fprintf('(%03d) theta0=%0.2f theta=%0.2f sigma0=%0.3f sigma=%0.3f\n',...
            k,theta0,theta,sigma0,sigma);

    % Run mean-field variational approximation.
    [logw(k),alpha(:,k),mu(:,k),s(:,k)] =...
    rss_varbvsr_bigmem_squarem(sumstats,sigbeta,logodds,options);

    clear sigbeta logodds theta0 theta sigma0 sigma;
    fprintf('\n');

  end

  clear data options gene_var;

end

