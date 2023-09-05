function out_mat = make_mat(in_lab,in_dat,in_par,in_opt,num_snp,num_gene)
% USAGE: create component matrices for a "network information" matrix
% INPUT:
%	in_lab: label of the component matrix, string
%	in_dat: data file of the component matrix, string
%	in_par: score mode of the component matrix, string
%	in_opt: score option of the component matrix, scalar
%	num_snp: total number of SNPs, scalar
%	num_gene: total number of genes, scalar
% OUTPUT:
%	out_mat: component matrix, sparse

  switch in_lab
    
    % Create SNP-gene proximity score matrix.
    case 'snp2gene'
      out_mat = make_snp2gene_mat(in_dat,in_par,in_opt,num_snp,num_gene);

    % Create gene-gene trans score matrix.
    case 'gene2gene'
      out_mat = make_gene2gene_mat(in_dat,in_par,in_opt,num_gene);

    % Create SNP-gene cis score matrix.
    case 'snp2gene_cis'
      out_mat = make_snp2gene_cis_mat(in_dat,in_par,in_opt,num_snp,num_gene);

    otherwise
      error('Invalid component matrix name ...');

  end

end

function snp2gene_mat = make_snp2gene_mat(snp2gene_file,snp2gene_par,snp2gene_opt,num_snp,num_gene)
% USAGE: output the sparse matrix of SNP-gene proximity scores
% INPUT:
%	snp2gene_file: a structure file with fields {rowid,colid,val,numsnp,numgene}, string
%	snp2gene_par: SNP-gene proximity score mode, string
%	snp2gene_opt: SNP-gene proximity score option, scalar
%	num_snp: total number of SNPs, scalar
%	num_gene: total number of genes, scalar
% OUTPUT:
%	snp2gene_mat: SNP by gene sparse matrix, SNP-gene proximity scores

  % Output an all zero sparse matrix when `snp2gene_file` is empty.
  if isempty(snp2gene_file)
    snp2gene_mat = sparse(num_snp,num_gene);
    return
  else
    snp2gene = matfile(snp2gene_file);
  end

  % Check the numbers of SNPs and genes.
  if (num_snp ~= snp2gene.numsnp)
    error('Inconsistent number of SNPs in `snp2gene` ...');
  end
  if (num_gene ~= snp2gene.numgene)
    error('Inconsistent number of genes in `snp2gene` ...');
  end

  % Confirm there is no duplicated SNP-gene pair.
  pos_tmp = [double(snp2gene.rowid) double(snp2gene.colid)];

  if (size(pos_tmp, 1) ~= size(unique(pos_tmp,'rows'), 1))
    error('At least one duplicated SNP-gene pair ...');
  end

  clear pos_tmp;

  % Create SNP-gene proximity scores.
  switch snp2gene_par

    % Set continuous SNP-gene proximity scores using an exponential model.
    % Apply an exponential decay to normalize SNP-gene distance to (0,1].
    % Use `snp2gene_opt` to control the decay rate (unit: base pair).
    case 'dist_exp'
      snp2gene_out = exp(-double(snp2gene.val) ./ snp2gene_opt);

    % Set binary (0/1) SNP-gene proximity scores using a genomic window.
    % Set 0 if the SNP-gene distance is greater than the given window size.
    % Use `snp2gene_opt` to control the window radius (unit: base pair).
    case 'dist_bin'
      snp2gene_out = double(double(snp2gene.val) <= snp2gene_opt);

    % Use pre-determined SNP-gene scores stored in `snp2gene.val`.
    case 'asis'
      if isscalar(snp2gene.val)
        snp2gene_out = double(snp2gene.val) * ones(length(snp2gene.rowid), 1);
      else
        snp2gene_out = double(snp2gene.val);
      end
 
    otherwise
      error('Invalid SNP-gene proximity score mode ...');

  end

  % Create the sparse matrix of SNP-gene proximity scores.
  snp2gene_mat = convert_mat(snp2gene_out,snp2gene.rowid,snp2gene.numsnp,...
                          snp2gene.colid,snp2gene.numgene);

  % Run sanity checks on SNP-gene proximity score matrix.
  check_size(snp2gene_mat,'snp2gene',snp2gene.numsnp,snp2gene.numgene);
  check_nnz(snp2gene_mat,'snp2gene',nonzeros(snp2gene_out));

  % Run additional sanity check on SNP-gene binary proximity scores.
  % See spones manual: https://www.mathworks.com/help/matlab/ref/spones.html.
  check_dist_bin = strcmp(snp2gene_par,'dist_bin');
  check_dist_bin = check_dist_bin && ~isequal(snp2gene_mat,spones(snp2gene_mat));
  if check_dist_bin
    error('Inconsistent SNP-gene binary proximity score ...');
  end

end

function gene2gene_mat = make_gene2gene_mat(gene2gene_file,gene2gene_par,gene2gene_opt,num_gene)
% USAGE: output the sparse matrix of gene-gene trans scores
% INPUT:
%	gene2gene_file: a structure file with fileds {rowid,colid,val,numgene}, string
%	gene2gene_par: gene-gene trans score mode, string
%	gene2gene_opt: gene-gene trans score option, scalar
%	num_gene: total number of genes, scalar
% OUTPUT:
%	gene2gene_mat: gene by gene sparse matrix, gene-gene trans scores

  % Output an all one diagonal matrix when `gene2gene_file` is empty.
  if isempty(gene2gene_file)
    gene2gene_mat = sparse(diag(ones(num_gene,1)));
    return
  else
    gene2gene = matfile(gene2gene_file);
  end

  % Check the numbers of genes.
  if (num_gene ~= gene2gene.numgene)
    error('Inconsistent number of genes in `gene2gene` ...');
  end

  % Confirm there is no duplicated gene-gene pair.
  pos_tmp = [double(gene2gene.rowid) double(gene2gene.colid)];

  if (size(pos_tmp, 1) ~= size(unique(pos_tmp,'rows'), 1))
    error('At least one duplicated gene-gene pair ...');
  end

  clear pos_tmp;

  % Ensure that input scores are between 0 and 1.
  check_val = all(double(gene2gene.val) >= 0);
  check_val = check_val & all(double(gene2gene.val) <= 1);

  if ~check_val
    error('Input gene-gene trans scores must be between 0 and 1 ...');
  end

  % Transform input gene-gene trans scores, if applicable.
  switch gene2gene_par

    % Use input gene-gene trans scores as is (default).
    case 'asis'
      gene2gene_out = double(gene2gene.val);

    % Set all gene-gene trans scores as zero (i.e. ignore trans effects).
    % This mode is mainly developed for "model assessment" simulations.
    case 'zero'

      disp('The mode `gene2gene_par=zero` is designed mainly for simulations ...');
      gene2gene_out = zeros(length(gene2gene.val), 1);
      gene2gene_out(double(gene2gene.rowid) == double(gene2gene.colid)) = 1;

    % Set all gene-gene trans scores as one (i.e. assume same and maximal trans effects).
    % This mode is mainly developed for "model assessment" simulations.
    case 'same'

      disp('The option `gene2gene_par=same` is designed mainly for simulations ...');
      gene2gene_out = ones(length(gene2gene.val), 1);

    % Randomly permute TF-TG pairs of an actual network (i.e. randomly connect true vertices).
    % The simulated and true networks have the same vertices but different edges.
    % This mode is mainly developed for "method comparison" simulations.
    case 'rand'

      disp('The mode `gene2gene_par=rand` is designed mainly for simulations ...');

      % Create a copy of the actual network.
      gene2gene_tmp = gene2gene;
      clear gene2gene;

      % Extract unique TFs and TGs from the actual network.
      gene2gene_colid = double(gene2gene_tmp.colid);
      gene2gene_rowid = double(gene2gene_tmp.rowid);

      tftg_flag = (gene2gene_rowid ~= gene2gene_colid);
      tftg_cols = unique(gene2gene_colid(tftg_flag));
      tftg_rows = unique(gene2gene_rowid(tftg_flag));

      % Create a grid of all possible TF-TG pairs.
      [X, Y]    = meshgrid(tftg_rows, tftg_cols);
      tftg_grid = [X(:) Y(:)];
      clear X Y;

      % Remove all ill-defined rows where TF==TG.
      tftg_grid = tftg_grid(tftg_grid(:,1) ~= tftg_grid(:,2), :);

      % Randomly draw unique TF-TG pairs from the grid without replacement.
      % Ensure that the number of TF-TG pairs matches the actual network.
      tftg_pair = datasample(1:size(tftg_grid,1),sum(tftg_flag),'Replace',false);

      % Ensure that all TFs and TGs appear in the simulated network,
      % because simulated and true networks must have the same vertices.
      check_cols = all(unique(tftg_grid(tftg_pair,2)) == unique(tftg_cols));
      check_rows = all(unique(tftg_grid(tftg_pair,1)) == unique(tftg_rows));

      % Repeat the permutation if not all TFs or TGs appear in previous permutation.
      while ~(check_cols && check_rows)
        tftg_pair = datasample(1:size(tftg_grid,1),sum(tftg_flag),'Replace',false);

        check_cols = all(unique(tftg_grid(tftg_pair,2)) == unique(tftg_cols));
        check_rows = all(unique(tftg_grid(tftg_pair,1)) == unique(tftg_rows));
      end
      clear tftg_cols tftg_rows check_cols check_rows; 

      % Replace the actual TF-TG pairs with permuted TF-TG pairs.
      gene2gene_colid(tftg_flag) = tftg_grid(tftg_pair, 2); 
      gene2gene_rowid(tftg_flag) = tftg_grid(tftg_pair, 1);
      clear tftg_flag tftg_grid tftg_pair;  

      % Create a TF-TG permuted version of the actual network.
      gene2gene.rowid   = gene2gene_rowid;
      gene2gene.colid   = gene2gene_colid;
      gene2gene.numgene = double(gene2gene_tmp.numgene);
      gene2gene_out     = double(gene2gene_tmp.val);
      clear gene2gene_tmp gene2gene_rowid gene2gene_colid;

    % Add random Gaussian noises to the input gene-gene trans scores.
    % This mode is mainly developed for "model perturbation" simulations.
    case 'add_noise'

      disp('The mode `gene2gene_par=add_noise` is designed mainly for simulations ...');
      gene2gene_out = double(gene2gene.val);

      % Add Gaussian noises to gene-gene trans scores.
      gene2gene_out = add_noise(gene2gene_out,gene2gene_opt);
 
      % Ensure gene-gene trans score matrix has correct diagonal values (1).
      gene2gene_out(double(gene2gene.rowid) == double(gene2gene.colid)) = 1;

    % Randomly set the input gene-gene trans scores to zero (i.e. remove TF-TG pairs).
    % This mode is mainly developed for "model perturbation" simulations.
    case 'cut_edge'

      disp('The mode `gene2gene_par=cut_edge` is designed mainly for simulations ...');
      gene2gene_out = double(gene2gene.val);

      % Find all TF-TG edges for random removal.
      tftg_flag = (double(gene2gene.rowid) ~= double(gene2gene.colid));

      % Randomly set gene-gene trans scores to zero, only for TF-TG pairs.
      gene2gene_out_tftg_only  = gene2gene_out(tftg_flag);
      gene2gene_out(tftg_flag) = cut_edge(gene2gene_out_tftg_only, gene2gene_opt);
      clear gene2gene_out_tftg_only;

      % Ensure gene-gene trans score matrix has correct diagonal values (1).
      if ~all(gene2gene_out(~tftg_flag) == 1)
        disp('Diagonal values of gene-gene trans score matrix must be all one ...');
      end

    otherwise
      error('Invalid gene-gene trans score mode ...');

  end

  % Create the sparse matrix of gene-gene trans scores.
  gene2gene_mat = convert_mat(gene2gene_out,gene2gene.rowid,gene2gene.numgene,...
                           gene2gene.colid,gene2gene.numgene);

  % Run sanity checks on gene-gene trans score matrix.
  check_size(gene2gene_mat,'gene2gene',gene2gene.numgene);
  check_nnz(gene2gene_mat,'gene2gene',nonzeros(gene2gene_out));
  check_diag(gene2gene_mat,'gene2gene',1);

end

function snp2gene_cis_mat = make_snp2gene_cis_mat(snp2gene_cis_file,snp2gene_cis_par,snp2gene_cis_opt,num_snp,num_gene)
% USAGE: output the sparse matrix of SNP-gene cis scores
% INPUT:
%	snp2gene_cis_file: a structure file with fileds {rowid,colid,val,numsnp,numgene}, string
%	snp2gene_cis_par: SNP-gene cis score mode, string
%	snp2gene_cis_opt: SNP-gene cis score option, scalar
%	num_snp: total number of SNPs, scalar
%	num_gene: total number of genes, scalar
% OUTPUT:
%	snp2gene_cis_mat: SNP by gene matrix, SNP-gene cis scores

  % Output an all zero sparse matrix when `snp2gene_cis_file` is empty.
  if isempty(snp2gene_cis_file)
    snp2gene_cis_mat = sparse(num_snp,num_gene);
    return
  else
    snp2gene_cis = matfile(snp2gene_cis_file);
  end

  % Check the numbers of SNPs and genes.
  if (num_snp ~= snp2gene_cis.numsnp)
    error('Inconsistent number of SNPs in `snp2gene_cis` ...');
  end
  if (num_gene ~= snp2gene_cis.numgene)
    error('Inconsistent number of genes in `snp2gene_cis` ...');
  end

  % Confirm there is no duplicated SNP-gene cis pair.
  pos_tmp = [double(snp2gene_cis.rowid) double(snp2gene_cis.colid)];

  if (size(pos_tmp, 1) ~= size(unique(pos_tmp,'rows'), 1))
    error('At least one duplicated SNP-gene cis pair ...');
  end

  clear pos_tmp;

  % Ensure input SNP-gene cis scores are between 0 and 1.
  check_val = all(double(snp2gene_cis.val) >= 0);
  check_val = check_val & all(double(snp2gene_cis.val) <= 1);

  if ~check_val
    error('Input SNP-gene cis scores must be between 0 and 1 ...');
  end

  % Transform input SNP-gene cis scores, if applicable.
  switch snp2gene_cis_par

    % Use input SNP-gene cis scores as is (default).
    case 'asis'
      snp2gene_cis_out = double(snp2gene_cis.val);

    % Set all SNP-gene cis scores as zero (i.e. ignore cis effects).
    % This mode is mainly developed for "model assessment" simulations.  
    case 'zero'

      disp('The mode `snp2gene_cis_par=zero` is designed mainly for simulations ...');
      snp2gene_cis_out = zeros(length(snp2gene_cis.val), 1);

    % Add random Gaussian noises to the input SNP-gene cis scores.
    % This mode is mainly developed for "model perturbation" simulations.
    case 'add_noise'

      disp('The mode `snp2gene_cis_par=add_noise` is designed mainly for simulations ...');
      snp2gene_cis_out = double(snp2gene_cis.val);

      % Add Gaussian noises to gene-gene trans scores.
      snp2gene_cis_out = add_noise(snp2gene_cis_out,snp2gene_cis_opt);

    % Randomly set the input SNP-gene cis scores to zero (i.e. remove RE-TG pairs).
    % This mode is mainly developed for "model perturbation" simulations.
    case 'cut_edge'

      disp('The mode `snp2gene_cis_par=cut_edge` is designed mainly for simulations ...');
      snp2gene_cis_out = double(snp2gene_cis.val);

      % Randomly set SNP-gene cis scores to zero.
      snp2gene_cis_out = cut_edge(snp2gene_cis_out,snp2gene_cis_opt);
    
    otherwise
      error('Invalid SNP-gene cis score mode ...');

  end

  % Create the sparse matrix of SNP-gene cis scores.
  snp2gene_cis_mat = convert_mat(snp2gene_cis_out,snp2gene_cis.rowid,snp2gene_cis.numsnp,...
                              snp2gene_cis.colid,snp2gene_cis.numgene);

  % Run sanity checks on SNP-gene cis score matrix.
  check_size(snp2gene_cis_mat,'snp2gene_cis',snp2gene_cis.numsnp,snp2gene_cis.numgene);
  check_nnz(snp2gene_cis_mat,'snp2gene_cis',nonzeros(snp2gene_cis_out));

end

function out_mat = convert_mat(out_val,row_id,row_num,col_id,col_num)
% USAGE: create a sparse matrix based on its nonzero entries

% NOTE: It is common to purposefully make the last line of the
% file include the desired size of the matrix with a value of 0.
% This practice ensures that the sparse matrix has desired size.
% See https://www.mathworks.com/help/matlab/ref/spconvert.html.

  % Create 3-column input matrix.
  out_long = zeros(length(out_val)+1, 3);

  out_long(:, 1) = [double(row_id); double(row_num)];
  out_long(:, 2) = [double(col_id); double(col_num)];
  out_long(:, 3) = [double(out_val); 0];

  % Confirm that there is no redundant row-column pair.
  rowcol = unique(out_long(1:end-1,1:2),'rows');
  if (size(rowcol,1) ~= length(out_val))
    error('At least one pair of row and column is redundant ...');
  end
  clear rowcol;

  % Convert 3-column input data to a sparse matrix.
  out_mat = spconvert(out_long);

end

function check_size(in_mat,in_name,row_num,col_num)
% USAGE: check the size for a given matrix

  if nargin == 3
    col_num = row_num;	  
  end

  row_num = double(row_num);
  col_num = double(col_num);

  if ~all(size(in_mat) == [row_num col_num])
    error(strcat(in_name,' matrix has incorrect size ...'));
  end

end

function check_nnz(in_mat,in_name,nz_val)
% USAGE: check the number of nonzero entries for a given matrix

  if (length(nz_val) ~= nnz(in_mat))
    error(strcat(in_name,' matrix has incorrect number of nonzero entries ...'));
  end	  

end

function check_diag(in_mat,in_name,diag_val)
% USAGE: check the diagonal values for a given matrix

  if ~all(diag(in_mat) == diag_val)
    error(strcat(in_name,' matrix has incorrect diagonal values ...'));
  end

end

function out_val = add_noise(in_val,in_opt)
% USAGE: add Gaussian noises to gene-gene trans or SNP-gene cis scores

  % Set standard deviation value for the Gaussian noise.
  % Use `in_opt` to control the relative size of noise.
  err_sd = in_opt * std(in_val);

  % Add Gaussian noises to the original scores.
  out_val = in_val + err_sd * randn(length(in_val), 1);

  % Ensure noise-added scores are still within [0,1] range.
  out_val = max(0, out_val);
  out_val = min(1, out_val);

end

function out_val = cut_edge(in_val,in_opt)
% USAGE: randomly set gene-gene trans or SNP-gene cis scores as zero

  out_val = in_val;

  % Sample trans gene-gene or cis SNP-gene pairs without replacement
  % such that all pairs have equal chances of being selected.
  % Use `in_opt` to control the proportion of removed pairs.
  pair_num = length(in_val);
  zero_num = round(in_opt * pair_num);
  zero_idx = datasample(1:pair_num,zero_num,'Replace',false);

  % Remove selected pairs by setting their scores as 0.
  out_val(zero_idx) = 0;

end
