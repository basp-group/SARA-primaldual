function [xsol, L1_v, L2_v, delta_v, sol_v, snr_v] = sdmm_bpconpar(y, epsilon, epsilons, A, At, T, W, Psi, Psit, param)
%
% [sol, L1_v, L2_v, delta_v, sol_v, snr_v] = sdmm_bpconpar(y, epsilon, epsilons, A, At, T, W, Psi, Psit, param)
% implements Algorithm 3 as described in [1].
% The implementation simulates a distributed setup for parallel processing.
%
% Inputs:
% y{:} - the visibility data
% epsilon - global L2 bound
% epsilons - global L2 stop criterion bound
% A - function handle, linear operator modeling the measurement process,
%     FFT and zero padding 
% At - the adjoint of A
% T{:} - gridding matrix containing the interpolation kernels (and any other
%     modeled DDEs)
% W{:} - selection of the FFT coefficients associated with each block of
%        data
% Psi - function handle, prior regularisation function
% Psit - the adjoint of Psi
% param - configuration parameters
% 
% Outputs: 
% xsol - final solution
% L1_v - evolution of the L1 norm
% L2_v - evolution of the L2 norm
% delta_v - evolution of the relative solution variation
% sol_v - solutions at inerations probvided in param.sol_steps
% snr_v - evolution of the SNR
%
% Authors: Rafael Carrillo & Alexandru Onose
%
% [1] A. Onose, R. E. Carrillo, A. Repetti, J. D. McEwen, J.-P. Thiran,
% J.-C. Pesquet and Y. Wiaux - Scalable splitting algorithms for big-data
% interferometric imaging in the SKA era, MNRAS (2016) 
% doi:10.1093/mnras/stw1859 

% number of nodes
R = length(y);
P = length(Psit);

% oversampling vectorized data length
No = size(W{1}, 1);
[Noy, Nox] = size(A(At(zeros(No, 1))));

% number of pixels
N = numel(At(zeros(No, 1)));
[Ny, Nx] = size(At(zeros(No, 1)));

% Optional input arguments
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'rel_obj'), param.rel_obj = 1e-4; end
if ~isfield(param, 'max_iter'), param.max_iter = 200; end
if ~isfield(param, 'gamma'), param.gamma = 1e-2; end
if ~isfield(param, 'max_iter_cg'), param.max_iter_cg = 20; end
if ~isfield(param, 'tol_cg'), param.tol_cg = 1e-6; end
if ~isfield(param, 'weights'), param.weights = 1; end
if ~isfield(param, 'sol_steps')
    param.sol_steps = inf;
else
    if param.sol_steps(end) ~= inf
        param.sol_steps = [param.sol_steps inf];
    end
end

% Useful functions for the projection
sc = @(z, bound) z*min(bound/norm(z(:)), 1); % scaling
soft = @(z, T) sign(z).*max(abs(z)-T, 0); %soft thresholding

% Initialization


if isfield(param, 'initsol')
    x3=param.initsol;
    ns = A(x3);
    for k = 1:R
        x2{k} =T{k}*ns(W{k});
        s2{k} =zeros(size(y{k}));
        z2{k} =zeros(size(y{k}));
    end
    x1=Psit(x3);
else
    x3=zeros(size(At(zeros(size(W{1})))));
    for k = 1:R
        x2{k} =zeros(size(y{k}));
        s2{k} =zeros(size(y{k}));
        z2{k} =zeros(size(y{k}));
    end
    x1=Psit(x3);  
end

z1=zeros(size(x1));
z3=zeros(size(x3));

% set up log variables
L1_v = zeros(param.max_iter, 1);
L2_v = zeros(param.max_iter, 1);

delta_v = zeros(param.max_iter, 1);
sol_steps = param.sol_steps;
sol_step_count = 1;
sol_v = zeros(length(sol_steps)-1, Ny, Nx);

snr_v = zeros(param.max_iter, 1);

%Linear operator for the cg solver
P = @(x) parallel_a(@(x)(1/sqrt(param.nu2)*A(x)),@(x) (1/sqrt(param.nu2)*At(x)),T,W,R,x);

%Initial solution estimate
xsol = x3;

flag = 0;
res1 = zeros(R,1);
epsilont = norm(epsilon);

%Sum reduce to compute data contribution
nst = zeros(No, 1);
nst(W{1}) = T{1}'*(x2{1} - z2{1});
g{1} = At(nst);
r = g{1};
for k = 2:R
    nst = zeros(No, 1);
    nst(W{k}) = T{k}'*(x2{k} - z2{k});
    g{k} = At(nst);
    r = r + g{k};
end

% Main loop
for t = 1:param.max_iter
    
    prev_xsol = xsol;
    norm_prevsol = norm(prev_xsol(:));
    
    %Mixing step. Final solution
    px = 1/param.nu1 * (Psi(x1-z1) + 1/param.nu2 * r + (x3-z3));
    xsol = solver_cg(P,px,0,param.max_iter_cg,xsol,param.tol_cg);
    
    
    
    % solution relative change
    if (norm_prevsol == 0)
        rel_sol_norm_change = 1;
    else
        rel_sol_norm_change = norm(xsol(:) - prev_xsol(:))/norm_prevsol;
    end
    
    s1 = Psit(xsol);
    curr_norm = sum(param.weights(:).*abs(s1(:)));
    
    %Proximal L1 operator
    %s1 = Psit(sol);
    x1 = soft(s1 + z1, param.gamma*param.weights);
    z1 = z1 + s1 - x1;
    
    ns = A(xsol);
    %Projection onto the L2-ball
    for k = 1:R
        s2{k} = T{k}*ns(W{k});
        res1(k) = norm(s2{k} - y{k});
        x2{k} = y{k} + sc(s2{k} + z2{k} - y{k}, epsilon(k));
        z2{k} = z2{k} + s2{k} - x2{k};
        nst = zeros(No, 1);
        nst(W{k}) = T{k}'*(x2{k} - z2{k});
        g{k} = At(nst);
    end
    
    %Projection on to the positive orthant
    x3 = real(xsol + z3);
    x3(x3<0) = 0;
    z3 = z3 + xsol - x3;
    
    %Log
    if (param.verbose >= 1)
        fprintf('Iter %i\n',t);
        fprintf(' L1 norm = %e, rel solution norm change = %e\n', ...
            curr_norm, rel_sol_norm_change);
        fprintf(' epsilon = %e, residual = %e\n\n', epsilont, norm(res1));
        if (param.verbose >= 2)
            for h = 1:R
                fprintf('  epsilon%i = %e, residual%i = %e\n\n', h, epsilon(h), h, res1(h));
            end
        end
    end
    if (param.verbose <= 0.5)
        fprintf('.\n');fprintf('\b');
        if mod(t, 50) == 0
            fprintf('\n');
        end
    end
    if (param.verbose >= 0.5)
        L1_v(t) = curr_norm;
        L2_v(t) = norm(res1);
        
        delta_v(t) = rel_sol_norm_change;
        snr_v(t) = 20*log10(norm(param.im(:))/norm(param.im(:) - xsol(:)));
        if t == sol_steps(sol_step_count)
            sol_v(sol_step_count, :, :) = xsol;
            sol_step_count = sol_step_count + 1;
        end
    end
    %Global stopping criteria
    if (rel_sol_norm_change < param.rel_obj && prod(res1 <= epsilons))
        flag = 1;
        break;
    end
    
    %Sum reduce to compute data contribution
    r = g{1};
    for k = 2:R
        r = r + g{k};
    end
    
end

%Final log
if (param.verbose > 0)
    if (flag == 1)
        fprintf('\n Solution found\n');
        fprintf(' Objective function = %e\n', curr_norm);
        fprintf(' Final residual = %e\n', res1);
    else
        fprintf('\n Maximum number of iterations reached\n');
        fprintf(' Objective function = %e\n', curr_norm);
        fprintf(' Relative variation = %e\n', rel_sol_norm_change);
        fprintf(' Final residual = %e\n', res1);
        fprintf(' epsilon = %e\n', epsilont);
    end
end

% trim the log vectors to the size of the actual iterations performed
if (param.verbose >= 0.5)
    L1_v = L1_v(1:t);
    L2_v = L2_v(1:t);
    delta_v = delta_v(1:t);
    sol_v = sol_v(1:sol_step_count-1, :, :);
    snr_v = snr_v(1:t);
end

end

function s = parallel_a(A,At,T,W,R,x)
    No = size(W{1}, 1);
    s = 2*x;
    ns = A(x);
    for k = 1:R
        nst = zeros(No, 1);
        nst(W{k}) = T{k}'*(T{k}*ns(W{k}));
        s = s + At(nst);
    end
end