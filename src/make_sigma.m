function [sigma0,sigma] = make_sigma(data,gene_var,logodds,eta,rho)
% USAGE: induce (sigma0, sigma) from (eta, rho) and data
% INPUT:
%	data: paths of input data files for RSS-NET, structure
%	gene_var: (net_mat.^2) * expr_nval, num_gene by 1
%	logodds: prior log odds ratio of each SNP, num_snp by 1
%	eta: a scalar between 0 and 1
%	rho: a scalar between 0 and 1
% OUTPUT:
%	sigma0: a non-negative scalar
%	sigma: a non-negative scalar

  % Extract the sample size of GWAS.
  num_sam = data.num_sam;

  % Load standard errors of single-SNP effect size estimates.
  gwas_data = matfile(data.sumstats_file);

  se  = cell2mat(gwas_data.se);
  ns2 = num_sam .* (se.^2);
  clear se;

  % Convert log odds to proportions.
  pival = 1 ./ (1 + exp(-logodds));

  % Compute scaling factors.
  bvec = pival ./ ns2;
  evec = (pival .* gene_var) ./ ns2;

  % Induce (sigma0, sigma) from (eta, rho).
  sigma0 = sqrt(eta*(1-rho) / sum(bvec));
  sigma  = sqrt(eta*rho / sum(evec));

end
