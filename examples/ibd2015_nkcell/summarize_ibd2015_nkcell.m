%%% PLEASE MODIFY PATHS BELOW IF NECESSARY %%%

src_path = '/home/users/xiangzhu/';
dat_path = '/scratch/users/xiangzhu/rss_net/ibd2015_nkcell/';
res_path = '/scratch/users/xiangzhu/rss_net/';

%%% PLEASE DO NOT MODIFY CODES BELOW IN GENERAL %%%

% specify regulatory network
net_name = 'Primary_Natural_Killer_cells_from_peripheral_blood'; 
cis_name = 'NK';

% specify complex trait
gwas_name = 'ibd2015';

% specify intermediate result path
res_path = strcat(res_path,gwas_name,'_',net_name,'_',cis_name,'/');

% specify hyper-parameter size
par_set = 125;

% run RSS-NET
run('summary_template.m');

