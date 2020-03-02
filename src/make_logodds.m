function logodds = make_logodds(snp2net_file,theta0,theta,options)
% USAGE: compute prior log odds based on SNP-network proximity
% INPUT:
%	snp2net_file: file of SNP-network proximity data, string
%	theta0 & theta: two RSS-NET hyper-parameters, scalar
%	options: user-specified behaviour of RSS-NET, structure
% OUTPUT:
%	logodds: log(prior PIP/(1-prior PIP)) for each SNP, num_snp by 1

  % Specify SNP-network proximity mode.
  if isfield(options,'snp2net_par')
    snp2net_par = options.snp2net_par;
  else
    snp2net_par = [];
  end

  % Load SNP-network proximity binary annotation vector.
  snp2net  = matfile(snp2net_file);
  net_flag = snp2net.val;
  clear snp2net;

  % Compute prior log odds based on SNP-network proximity.
  % Note that (theta0, theta) are on log base 10,
  % and `logodds` in RSS programs uses natural logarithm.
  % See https://github.com/stephenslab/rss/blob/master/src_vb/rss_varbvsr_bigmem_squarem.m.
  logodds = zeros(length(net_flag), 1);

  logodds(net_flag==0) = log(10) * theta0;
  logodds(net_flag==1) = log(10) * (theta0+theta);

  % Ignore SNP-network proximity binary annotation.
  % This mode is mainly developed for "model assessment" simulations.
  if strcmp(snp2net_par,'zero')
    logodds = (log(10) * theta0) * ones(length(net_flag), 1);
    disp('Ignore SNP-network proximity binary annotation ...');
  end

end
