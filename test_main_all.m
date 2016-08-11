close all;

addpath data/
addpath lib/
addpath alg/
try
    run('../../../lib/irt/setup.m');
catch
    fprintf('IRT lib not found \nSet the path to irt/setup.m in test_main_all.m\n');
    break;
end



%% run parameters
% 0 - loads new data from file based on the dataset number supplied
% 1 - generates new data
% 2 - uses the data in matlab's workspace
gen_data = 1;
gen_figures = 1;
gen_only_average_figures = 1;


save_dataset_number = 0; % number of the dataset to write files to
save_dataset_subnumber = 0; % number of the dataset to write files to

save_data_on_disk = 1; % flag
save_path = 'results/';

% number of tests to perform on the same coverage with different noise
% realisations
num_tests = 1;


run_pdfb_bpcon_par_sim_rescaled = 1; % flag
run_pdfb_bpcon_par_sim_rand_rescaled = 1; % flag
run_admm_bpconpar = 1; %flag
run_sdmm_bpconpar = 1; %flag


%% input test file name
file_name = 'M31.fits';
% file_name = 'W28.fits';
% file_name = 'CYN.fits';

%% config parameters
verbosity = 1;

input_snr = 20; % noise level (on the measurements)

nlevel = 4; % wavelet level
ox = 2; % oversampling factors for nufft
oy = 2; % oversampling factors for nufft
Kx = 8; % number of neighbours for nufft
Ky = 8; % number of neighbours for nufft

%             L2 ball definition    stop criterion
eps_choice = [2                     3; % number of sigmas' above the mean for the L2 ball and the stopping criterion respectivelly
              0.01                  0.001; % p_val - 1-CDF the L2 ball and the stopping criterion respectivelly
              0                     1.005; % percentage of the ball epsilon
              1                     1]; % choice between sigmas if 1, p_val if 2 and percentage if 3
wlt_basis = {'db1', 'db2', 'db3', 'db4', 'db5', 'db6', 'db7', 'db8', 'self'}; % wavelet basis to be used
use_same_stop_criterion = 1; % forces the distributed criterion to have same norm as in the nondistributed setup
                 
%% samplig pattern parameters
% options 'gaussian', 'file', 'gaussian+large-holes', 'file+undersample'

% sampling_pattern = 'file+undersample';
sampling_pattern = 'gaussian';
% sampling_pattern = 'gaussian+large-holes';
% sampling_pattern = 'file';

sparam.file_name = 'data/test_vla.uvw.mat'; % file name for uv coverage
sparam.p = 1; % number of measurements as proportion of number of pixels to recover

% for randomly generated coverages
sparam.hole_number = 8000; % number of holes to introduce for 'gaussian+large-holes'
sparam.hole_size = pi/60; % size of the missing frequency data
sparam.sigma = pi/4; % variance of the gaussion over continous frequency
sparam.sigma_holes = pi/3; % variance of the gaussion for the holes

% data splitting config
sparam.equal_partitioning = 1; % flag
sparam.equal_partitioning_no = 4;

% if sparam.equal_partitioning = 0 the partitioning is given by sparam.fpartition
sparam.fpartition = [icdf('norm', 0.25, 0, sparam.sigma), 0, icdf('norm', 0.75, 0, sparam.sigma), pi]; % partition (symetrically) the data to nodes (frequency ranges)


%% get input data
script_get_input_data;


%% PDFB parameter structure sent to the algorithm
param_pdfb.im = im; % original image, used to compute the SNR
param_pdfb.verbose = verbosity; % print log or not
param_pdfb.nu1 = 1; % bound on the norm of the operator Psi
param_pdfb.nu2 = evl; % bound on the norm of the operator A*G
param_pdfb.gamma = 1e-3; % convergence parameter L1 (soft th parameter)
param_pdfb.tau = 0.49; % forward descent step size
param_pdfb.rel_obj = 0; % stopping criterion
param_pdfb.max_iter = 50; % max number of iterations
param_pdfb.sol_steps = [inf]; % saves images at the given iterations


%% PDFB PROB parameter structure sent to the algorithm 
param_pdfb_prob.im = im; % original image, used to compute the SNR
param_pdfb_prob.verbose = verbosity; % print log or not
param_pdfb_prob.nu1 = 1; % bound on the norm of the operator Psi
param_pdfb_prob.nu2 = evl; % bound on the norm of the operator A*G
param_pdfb_prob.gamma = 1e-3; % convergence parameter L1 (soft th parameter)
param_pdfb_prob.tau = 0.49; % forward descent step size
param_pdfb_prob.rel_obj = 0; % stopping criterion
param_pdfb_prob.max_iter = 50; % max number of iterations
param_pdfb_prob.pd1 = 1; % probability of the selection of each basis at each iteration
param_pdfb_prob.pd2 = 0.75; % probability of the selection of each input data split at each iteration

%% ADMM parameter structure sent to the algorithm
param_admm.im = im; % original image, used to compute the SNR
param_admm.verbose = verbosity; % print log or not
param_admm.gamma = 1e-3*evl; % convergence parametter
param_admm.nu = evl; % bound on the norm of the operator A*G
param_admm.rel_obj = 0; % stopping criterion
param_admm.max_iter = 50; % max number of iterations
param_admm.sol_steps = [inf]; % saves images at the given iterations

% L1 sub-iterations parameters
param_admm.tight_L1 = 0; % indicate if Psit is a tight frame (1) or not (0)
param_admm.max_iter_L1 = 100; % max number of iterations to be performed for estimating the L1 prox
param_admm.rel_obj_L1 = 1e-3; % stopping criterion for the L1 prox
param_admm.pos_L1 = 1; % constrain for the positivity
param_admm.nu_L1 = 1; % bound on the norm of the operator Psi
param_admm.verbose_L1 = 0; % print log or not



%% SDMM parameter structure sent to the algorithm
param_sdmm.im = im; % original image, used to compute the SNR
param_sdmm.verbose = verbosity; % print log or not
param_sdmm.gamma = 1e-3; % convergence parameter
param_sdmm.rel_obj = 0; % stopping criterion
param_sdmm.max_iter = 50; % max number of iterations
param_sdmm.max_iter_cg = 50; % max number of iterations for conjugate gradient
param_sdmm.tol_cg = 1e-3; % conjugate gradient stopping criterion
param_sdmm.nu1 = 1; % bound on the norm of the operator Psi
param_sdmm.nu2 = evl; % bound on the norm of the operator A*G
param_sdmm.sol_steps = [inf]; % saves images at the given iterations


%% compute the solution
fprintf('Starting algorithms:\n\n');
tstart = tic;

script_run_all_tests_serial;

tend = toc(tstart);
fprintf('All algorithms runtime: %ds\n\n', ceil(tend));
