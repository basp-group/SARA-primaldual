% script that loads or generates the input data for the test
% if save_data_on_disk = 1 it saves the data on disk 


if gen_data == 0
    fprintf('Loading data from disk ... \n\n');
    if exist(sprintf('%s%s_input_data.mat', save_path, int2str(save_dataset_number)), 'file') || ... 
          exist(sprintf('%s%s_input_data_config.mat', save_path, int2str(save_dataset_number)), 'file')
        load(sprintf('%s%s_input_data.mat', save_path, int2str(save_dataset_number)));
        load(sprintf('%s%s_input_data_config.mat', save_path, int2str(save_dataset_number)));

        R = length(v);
        
        [A, At, G, W, Gw] = op_p_nufft([v u], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2]);
        T = G;
        Tw = Gw;
        
        [Psi, Psit] = op_p_sp_wlt_basis(wlt_basis, nlevel, Ny, Nx);
        [Psiw, Psitw] = op_sp_wlt_basis(wlt_basis, nlevel, Ny, Nx);
        
        %% compute the operator norm
        
        fprintf('Computing operator norm ...\n');
        evl = op_norm(@(x) Gw * A(x), @(x) At(Gw' * x), [Ny, Nx], 1e-4, 200, 2);

        %% L2 ball sizes
        epsilon = cell(num_tests, 1);
        epsilons = cell(num_tests, 1);
        epsilonT = cell(num_tests, 1);
        epsilonTs = cell(num_tests, 1);
        %% generate bounds
        for k = 1:num_tests
            [epsilonT{k}, epsilonTs{k}, epsilon{k}, epsilons{k}] = util_gen_L2_bounds(y0{k}, input_snr, eps_choice, use_same_stop_criterion);
        end

    else
        error('Unable to find data files');
    end
end

if gen_data == 1
    fprintf('Generating new data ... \n\n');
    %% image loading
    [im, N, Ny, Nx] = util_read_image(file_name);

    %% generate the sampling pattern
    sparam.N = N; % number of pixels in the image
    sparam.Nox = ox*Nx; % number of pixels in the image
    sparam.Noy = oy*Ny; % number of pixels in the image
    [u, v, ~, ~, ~, uvidx] = util_gen_sampling_pattern(sampling_pattern, sparam);


    %% measurement operator initialization 
    fprintf('Initializing the NUFFT operator\n\n');
    tstart = tic;
    
    [A, At, G, W, Gw] = op_p_nufft([v u], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2]);
    T = G;
    Tw = Gw;
    tend = toc(tstart);
    
    fprintf('Initialization runtime: %ds\n\n', ceil(tend));
    R = length(v);
    
    y0 = cell(num_tests, 1);
    y = cell(num_tests, 1);
    aY = cell(num_tests, 1);
    epsilon = cell(num_tests, 1);
    epsilons = cell(num_tests, 1);
    epsilonT = cell(num_tests, 1);
    epsilonTs = cell(num_tests, 1);
    y0f = cell(num_tests, 1);
    yf = cell(num_tests, 1);
    
    %% generate noisy input data
    for k = 1:num_tests
        [y0{k}, y{k}, y0f{k}, yf{k}, aY{k}] = util_gen_input_data(im, G, W, A, input_snr);
        [epsilonT{k}, epsilonTs{k}, epsilon{k}, epsilons{k}] = util_gen_L2_bounds(y0{k}, input_snr, eps_choice, use_same_stop_criterion);
    end
    yT = y;

    %% sparsity operator definition

    [Psi, Psit] = op_p_sp_wlt_basis(wlt_basis, nlevel, Ny, Nx);
    [Psiw, Psitw] = op_sp_wlt_basis(wlt_basis, nlevel, Ny, Nx);
    
    
    %% compute the operator norm
    fprintf('Computing operator norm ...\n');
    evl = op_norm(@(x) Gw * A(x), @(x) At(Gw' * x), [Ny, Nx], 1e-4, 200, verbosity);

    
    %% save data
    if save_data_on_disk == 1
        fprintf('Saving new data ... \n');

        if save_data_on_disk
            file_start = save_dataset_number;
            while exist(sprintf('%s%s_input_data.mat', save_path, int2str(file_start)), 'file') || ... 
                  exist(sprintf('$s%s_input_data_config.mat', save_path, int2str(file_start)), 'file')
                file_start = file_start + 1;
            end
            if file_start ~= save_dataset_number;
                fprintf('WARNING: Saving new data in file %d instead of %d \n\n', file_start, save_dataset_number);
            end

            save(sprintf('%s%s_input_data', save_path, int2str(file_start)), '-v7.3', 'im', 'y0', 'y', ... % 'G', 'Gw', 'W'
                'y0f', 'yf', 'evl'); % do not save G, it is large and is faster to compute than saving 
            save(sprintf('%s%s_input_data_config', save_path, int2str(file_start)), 'nlevel', 'N', 'Ny', ...
                'Nx', 'Ky', 'Kx', 'oy', 'ox', 'u', 'v', 'sampling_pattern', 'sparam', 'input_snr', 'file_name', 'wlt_basis', 'num_tests');
            
            if strcmp(sampling_pattern, 'file')
                for k = 1:num_tests
                    y_ = yf{k};
                    y0_ = y0f{k};
                    ccc = [uvidx{:}];
                    vis = zeros(size(y_));
                    vis(ccc(:)) = y_;
                    save(sprintf('%s%s_vis_data_t%s', save_path, int2str(file_start), int2str(k)), 'vis');
                    vis = zeros(size(y_));
                    vis(ccc(:)) = y_ - y0_;
                    save(sprintf('%s%s_vis_data_noise_t%s', save_path, int2str(file_start), int2str(k)), 'vis');
                    clear vis y_ y0_;
                end
            end
        end
        
    end
end

if gen_data == 2
    fprintf('Using data from workspace ... \n\n');
end
