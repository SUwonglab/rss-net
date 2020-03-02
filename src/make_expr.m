function expr_nval = make_expr(data,options)
% USAGE: create a vector of normalized gene expression values
% INPUT:
%	data: paths of input data files for RSS-NET, structure
%	options: user-specified behaviour of RSS-NET, structure
% OUTPUT:
%	expr_nval: normalized gene expression values, num_gene by 1

  % Output an all one vector if the expression is empty.
  % Currently this is the most recommended option.
  expr_file = data.expression_file;

  if isempty(expr_file)
    expr_nval = ones(data.num_gene, 1);
    return
  else
    warning('Currently RSS-NET does not know how to exploit gene expression effectively.');
    disp('Please either i) input an empty expression vector or ii) proceed with extreme caution.');
  end

  % NOTE:
  % Due to the `return` function above, the following codes are rarely used.
  % Please read the warning messages above.
  % I am still working on how to effectively exploit gene expression in RSS-NET.
 
  % Load expression-related mode.
  if isfield(options,'expression_par')
    expr_par = options.expression_par;
  else
    expr_par = [];
  end

  % Load gene expression data.
  expr = matfile(expr_file);

  % Confirm there is no duplicated gene.
  if (length(expr.genenid) ~= length(unique(expr.genenid)))
    error('At least one duplicated gene ...\n');
  end

  % Obtain z-transformed expression values.
  expr_val = expr.z;

  % Normalize z-transformed expression values.
  expr_nval = 1 + max(0, expr_val);

  % Force all genes to have the same expression level.
  % This mode is mainly developed for "model assessment" simulations.
  if strcmp(expr_par,'same')
    expr_nval = ones(size(expr_val));
    disp('Ignore context-specific expression pattern ...');
  end

end
