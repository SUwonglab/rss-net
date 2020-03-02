function net_mat = make_netmat(data,options)
% USAGE: create "network information" matrix for a given regulatory network
% INPUT:
%	data: input data files to create "network information" matrix, structure
%		- snp2gene_file: a mat file with {rowid,colid,val,numsnp,numgene}, string
%		- gene2gene_file: a mat file with {rowid,colid,val,numgene}, string
%		- snp2gene_cis_file: a mat file with {rowid,colid,val,numsnp,numgene}, string
%		- num_snp: total number of SNPs, scalar
%		- num_gene: total number of genes, scalar
%	options: user-specified parameters, structure 
%		- snp2gene_par: SNP-gene score mode, string
%		- gene2gene_par: gene-gene score mode, string
%		- snp2gene_cis_par: SNP-gene cis score mode, string
%		- snp2gene_opt: SNP-gene score option (mainly for simulations), scalar
%		- gene2gene_opt: gene-gene score option (mainly for simulations), scalar
%		- snp2gene_cis_opt: SNP-gene cis score option (mainly for simulations), scalar
% OUTPUT:
%	net_mat: "network information" sparse matrix, num_snp by num_gene

  % Check input arguments.
  if ~exist('data','var')
    error('Input data files are missing ...');
  end

  if ~exist('options','var')
    options = [];
  end

  % Specify input data file names.

  % row: SNP; col: gene; val: SNP-gene proximity score
  % source: snp2gene.R
  if isfield(data,'snp2gene_file')
    snp2gene_file = data.snp2gene_file;
  else
    error('snp2gene_file is missing ...');
  end

  % row: gene; col: gene; val: gene-gene trans score
  % source: gene2gene.R
  if isfield(data,'gene2gene_file')
    gene2gene_file = data.gene2gene_file;
  else
    error('gene2gene_file is missing ...');
  end

  % row: SNP; col: gene; val: SNP-gene cis score
  % source: snp2gene_cis.R
  if isfield(data,'snp2gene_cis_file')
    snp2gene_cis_file = data.snp2gene_cis_file;
  else
    error('snp2gene_cis_file is missing ...');
  end

  % Load and display important user-specified parameters.

  if isfield(options,'snp2gene_par')
    snp2gene_par = options.snp2gene_par;
  else
    snp2gene_par = 'dist_bin';
    snp2gene_opt = 1e6;
  end
  fprintf('SNP-gene proximity score mode: %s ...\n', snp2gene_par);

  if isfield(options,'gene2gene_par')
    gene2gene_par = options.gene2gene_par;
  else
    gene2gene_par = 'asis';
    gene2gene_opt = [];
  end
  fprintf('Gene-gene trans score mode: %s ...\n', gene2gene_par);

  if isfield(options,'snp2gene_cis_par')
    snp2gene_cis_par = options.snp2gene_cis_par;
  else
    snp2gene_cis_par = 'asis';
    snp2gene_cis_opt = [];
  end
  fprintf('SNP-gene cis score mode: %s ...\n', snp2gene_cis_par);

  % Load additional parameters developed mainly for simulations.

  if isfield(options,'snp2gene_opt')
    snp2gene_opt = options.snp2gene_opt;
  end

  if isfield(options,'gene2gene_opt')
    gene2gene_opt = options.gene2gene_opt;
  end

  if isfield(options,'snp2gene_cis_opt')
    snp2gene_cis_opt = options.snp2gene_cis_opt;
  end

  % Get the total numbers of SNPs and genes.
  num_snp  = data.num_snp;
  num_gene = data.num_gene;

  % Construct the sparse matrix of SNP-gene proximity scores.
  tic;
  snp2gene_mat = make_mat('snp2gene',snp2gene_file,snp2gene_par,snp2gene_opt,num_snp,num_gene);
  snp2gene_time=toc;

  fprintf('Time to create the SNP-gene proximity score matrix: %8.3f seconds ...\n', snp2gene_time);
  clear snp2gene_file snp2gene_time snp2gene_par snp2gene_opt;

  % Construct the sparse matrix of gene-gene trans scores.
  tic;
  gene2gene_mat = make_mat('gene2gene',gene2gene_file,gene2gene_par,gene2gene_opt,num_snp,num_gene);
  gene2gene_time=toc;

  fprintf('Time to create the gene-gene trans score matrix: %8.3f seconds ...\n', gene2gene_time);
  clear gene2gene_file gene2gene_time gene2gene_par gene2gene_opt;

  % Construct the sparse matrix of SNP-gene cis scores.
  tic;
  snp2gene_cis_mat = make_mat('snp2gene_cis',snp2gene_cis_file,snp2gene_cis_par,snp2gene_cis_opt,num_snp,num_gene);
  snp2gene_cis_time=toc;

  fprintf('Time to create the SNP-gene cis score matrix: %8.3f seconds ...\n', snp2gene_cis_time);
  clear snp2gene_cis_file snp2gene_cis_time snp2gene_cis_par snp2gene_cis_opt;

  % Compute the "network information" sparse matrix.
  tic;

  % Construct the total sparse matrix of SNP-gene scores.
  % total SNP-gene score := SNP-gene proximity score * (1 + SNP-gene cis score)
  snp2gene_all_mat = snp2gene_mat .* (spones(snp2gene_mat) + snp2gene_cis_mat);
  clear snp2gene_cis_mat snp2gene_mat;

  % SNP-gene network matrix := SNP-gene matrix %*% gene-gene matrix
  net_mat = snp2gene_all_mat * gene2gene_mat;
  nm_time = toc;

  fprintf('Time to compute the network information matrix: %8.3f seconds ...\n', nm_time);
  clear snp2gene_all_mat gene2gene_mat nm_time;

  % Double check the dimensions of network sparse matrix.
  if (size(net_mat,1) ~= num_snp)
    error('Inconsistent row numbers (num_snp) of network matrix ...');
  end
  if (size(net_mat,2) ~= num_gene)
    error('Inconsistent column numbers (num_gene) of network matrix ...');
  end

end

