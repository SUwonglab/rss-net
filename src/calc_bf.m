function [log10_bf,log10_bf_nt,log10_bf_ns,log10_bf_ts] = calc_bf(logw,param)
% USAGE: compute enrichment Bayes factors (BFs) using RSS-NET model fitting results
% INPUT:
%	logw: log variational lower bound (up to a constant), num_hyper by 1
%	param: hyper-parameter grid [theta0 theta sigma0 sigma], num_hyper by 4
% OUTPUT:
%	log10_bf: log10 BF of the full model with theta > 0 or sigma > 0, scalar
%	log10_bf_nt: log10 BF of the reduced model with theta == 0, scalar
%	log10_bf_ns: log10 BF of the reduced model with sigma == 0, scalar
%	log10_bf_ts: log10 BF of the reduced model with theta > 0 and sigma > 0, scalar

% NOTE:
% The log10 BF of the reduced model with sigma == 0, `log10_bf_ns`,
% is the same as the log10 enrichment BF based on RSS-E model.
% See https://doi.org/10.1038/s41467-018-06805-x for more details. 

  % Define the reduced model with theta == 0.
  nt = (param(:,2) == 0);

  % Define the reduced model with sigma == 0.
  ns = (param(:,4) == 0);

  % Define the reduced model with theta > 0 and sigma > 0.
  ts = ((~nt) & (~ns));

  % Define the baseline model with theta == 0 and sigma == 0.
  nn = (nt & ns);

  % Estimate log marginal likelihood of the baseline model.
  loglik_nn = avg_logw(logw(nn));

  % Estimate log marginal likelihoods of three reduced models.
  loglik_nt = avg_logw(logw(nt & (~ns)));
  loglik_ns = avg_logw(logw(ns & (~nt)));
  loglik_ts = avg_logw(logw(ts));

  % Estimate log marginal likelihood of the full model.
  loglik = avg_logw(logw);

  % Estimate log10 BFs under one full and three reduced models.
  log10_bf    = (loglik-loglik_nn) / log(10);
  log10_bf_nt = (loglik_nt-loglik_nn) / log(10);
  log10_bf_ns = (loglik_ns-loglik_nn) / log(10);
  log10_bf_ts = (loglik_ts-loglik_nn) / log(10);

end

function loglik = avg_logw(logw)
% USAGE: estimate the marginal log-likelihood from log variational lower bound
% INPUT:
%	logw: log variational lower bound (up to a constant), num_hyper by 1
% OUTPUT:
%	loglik: estimated marginal log-likelihood, scalar
% SOURCE: https://github.com/pcarbo/varbvs/blob/master/varbvs-MATLAB/bayesfactor.m

  % Find the largest entry of logw as a common constant.
  c = max(logw(:));

  % Average all entries of logw to compute marginal log-likelihood.
  loglik = c + log(mean(exp(logw(:) - c)));

end
