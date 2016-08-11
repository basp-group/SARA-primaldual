function [y0, y, y0f, yf, aY] = util_gen_input_data(im, T, W, A, input_snr)
% generates the input data


R = length(T);

y0f = A(im);
y0 = cell(R, 1);
for q = 1:R
    y0{q} = T{q} * y0f(W{q});
end

Nm = numel(cell2mat(y0));
normy0 = norm(cell2mat(y0));

% add Gaussian i.i.d. noise
sigma_noise = 10^(-input_snr/20) * normy0/sqrt(Nm);

y = cell(R, 1);
aY = cell(R, 1);
for q = 1:R
    Nmq = length(y0{q});
    noise = (randn(Nmq, 1) + 1i*randn(Nmq, 1)) * sigma_noise/sqrt(2);
    y{q} = y0{q};
    y{q} = y{q} + noise;
    aY{q} = abs(y{q})./sigma_noise;
end
yf = cell2mat(y);
y0f = cell2mat(y0);

end

