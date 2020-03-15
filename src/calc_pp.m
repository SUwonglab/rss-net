function [pp,pp_nt,pp_ns,pp_nn] = calc_pp(seg_path,logw,alpha,param)
% USAGE: compute locus-level association statistics using RSS-NET model fitting results
% INPUT:
%	seg_path: mat file path of genomic locus definition, string
%	logw: log variational lower bound (up to a constant), num_hyper by 1
%	alpha: variational posterior inclusion probabilities, num_snp by num_hyper
%	param: hyper-parameter grid [theta0 theta sigma0 sigma], num_hyper by 4
% OUTPUT:
%	pp: P1 and P2 under the full model, num_locus by 1
%	pp_nt: P1 and P2 under the reduced model with theta == 0, num_locus by 1
%	pp_ns: P1 and P2 under the reduced model with sigma == 0, num_locus by 1
%	pp_nn: P1 and P2 under the baseline model with theta == 0 and sigma == 0, num_locus by 1

% NOTE:
% The locus-level P1 and P2 are defined as the posterior probabilities that
% at least one SNP and two SNPs in a locus is trait-associated, respectively.
% See https://doi.org/10.1038/s41467-018-06805-x for more details.
% Calculation of P1 and P2 statistics is implemented as follows:
% https://github.com/stephenslab/rss/blob/master/src_vb/compute_pip.m. 

  % Preallocate locus-level P1 and P2 output.
  seg_info = matfile(seg_path);
  num_segs = size(seg_info.Aseg, 2);

  pp    = zeros(num_segs,2);
  pp_nt = zeros(num_segs,2);
  pp_ns = zeros(num_segs,2);
  pp_nn = zeros(num_segs,2);

  clear seg_info num_segs;

  % Define the reduced model with theta == 0.
  nt = (param(:,2) == 0);

  % Define the reduced model with sigma == 0.
  ns = (param(:,4) == 0);

  % Define the baseline model with theta == 0 and sigma == 0.
  nn = (nt & ns);

  % Estimate locus-level P1 and P2 under the baseline model.
  [pp_nn(:,1),pp_nn(:,2)] = compute_pip(seg_path,logw(nn),alpha(:,nn));

  % Estimate locus-level P1 and P2 under two reduced models.
  [pp_nt(:,1),pp_nt(:,2)] = compute_pip(seg_path,logw(nt),alpha(:,nt));
  [pp_ns(:,1),pp_ns(:,2)] = compute_pip(seg_path,logw(ns),alpha(:,ns));
  
  % Estimate locus-level P1 and P2 under the full model.
  [pp(:,1),pp(:,2)] = compute_pip(seg_path,logw,alpha);

end
