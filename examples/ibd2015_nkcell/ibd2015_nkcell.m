%%% PLEASE MODIFY PATHS BELOW IF NECESSARY %%%

src_path = '/home/users/xiangzhu/';
dat_path = '/scratch/users/xiangzhu/rss_net/ibd2015_nkcell/';
out_path = '/scratch/users/xiangzhu/rss_net/';

%%% PLEASE DO NOT MODIFY CODES BELOW IN GENERAL %%%

% specify regulatory network
net_name = 'Primary_Natural_Killer_cells_from_peripheral_blood'; 
cis_name = 'NK';

% specify total number of autosome protein-coding genes
ngene = 18334;

% specify complex trait
gwas_name = 'ibd2015';

% specify GWAS sample size
nsam = (12882+21770);

% specify hyper-parameters
eta_set    = 0.3;
rho_set    = (0:0.2:0.8);
theta0_set = (-3:0.05:-2.8);
theta_set  = (0:0.25:1);

hyper_data.theta0 = theta0_set;
hyper_data.theta  = theta_set;
hyper_data.eta    = eta_set;
hyper_data.rho    = rho_set;

clear *_set;

% run RSS-NET
run('analysis_template.m');

% exit the program
exit;
