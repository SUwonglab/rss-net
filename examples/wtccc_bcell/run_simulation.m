%%% PLEASE MODIFY PATHS BELOW IF NECESSARY %%%

src_path = '/home/xiangzhu/';
dat_path = '/scratch/midway2/xiangzhu/tmp/rss_net/wtccc_net/';
out_path = '/scratch/midway2/xiangzhu/tmp/rss_net/wtccc_net_out/';

%%% PLEASE DO NOT MODIFY CODES BELOW IN GENERAL %%%

true_theta0 = -4;
true_theta  = 2;
true_sigb   = 1;
true_pve    = 0.6;

switch trial_name
  case 'm0'
    true_sige = 0;
  case 'm1'
    true_sige = 2;
end

% specify random seed
seed = 2104;

% test the simulation template
run('simulation_template.m');

% exit the program
exit;

