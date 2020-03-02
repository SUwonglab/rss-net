function [hyper_param,induced_flag] = make_hpgrid(hyper)
% USAGE: create a grid of all possible hyper-parameter combinations for RSS-NET
% INPUT:
%	hyper: hyper-parameter values for RSS-NET, structure
% OUTPUT:
%	hyper_param: a grid of all possible hyper-parameter combinations, num_hyper by 4
%	induced_flag: true if prior of (sigma0, sigma) is induced from prior of (eta, rho), logical 

  % Load the set of `theta0` values.
  if isfield(hyper,'theta0')
    theta0_set = hyper.theta0;
  else
    error('Hyper-parameter `theta0` must be provided ...');
  end

  % Load the set of `theta` values.
  if isfield(hyper,'theta')
    theta_set = hyper.theta;
  else
    error('Hyper-parameter `theta` must be provided ...');
  end

  % Determine whether to induce prior for effect size.
  % a) param == [theta0 theta eta rho] if induced_flag == true
  % b) param == [theta0 theta sigma0 sigma] if induced_flag == false

  eta_n_rho = isfield(hyper,'eta') & isfield(hyper,'rho');
  if eta_n_rho
    induced_flag = true; %#ok<*NASGU>

    third_set = hyper.eta;
    forth_set = hyper.rho;
  end

  sigma0_n_sigma = isfield(hyper,'sigma0') & isfield(hyper,'sigma');
  if sigma0_n_sigma
    induced_flag = false;

    third_set = hyper.sigma0;
    forth_set = hyper.sigma;
  end

  if ~(eta_n_rho || sigma0_n_sigma)
    error('Hyper-parameters `[eta,rho]` or `[sigma0,sigma]` must be provided ...');
  end

  if (eta_n_rho && sigma0_n_sigma)
    error('Hyper-parameters `[eta,rho]` and `[sigma0,sigma]` cannot coexist ...');
  end

  % Create a 4-D rectangular grid.
  % See `ndgrid` manual: https://www.mathworks.com/help/matlab/ref/ndgrid.html.
  [X1,X2,X3,X4] = ndgrid(theta0_set,theta_set,third_set,forth_set);

  % Reshape 4-D to 2-D grid (i.e. 4-column matrix).
  hyper_param = [X1(:) X2(:) X3(:) X4(:)];
  clear X1 X2 X3 X4;

end
