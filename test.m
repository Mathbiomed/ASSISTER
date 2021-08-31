prop_functions = @propensity_functions;
Stoi = [-1 -1 1 0 0 ; 1 1 -1 0 0 ; 0 0 0 1 0;...
         0 0 0 0 1; 0 0 0 -1 0; 0 0 0 0 -1]';
tspan = 0:200:1600;
err_tol = 0.1;
x_init = [10,10,0,0,0]';
rev_idx = [1 2];

species = Gillespie_Reduction(Stoi, prop_functions, x_init, tspan, rev_idx, err_tol);

function lambda = propensity_functions(x)
%fast reaction const
k_f = 100; %binding
k_b = 1; %unbinding
%slow reaction const
k_R = 0.1; %M_R transcription
k_A = 0.001; %M_A transcription
k_Rd = 0.001; %M_R Degradation
k_Ad = 0.001; %M_A Degradation
lambda(1) = k_f * x(1) * x(2); % D + P -> D:P
lambda(2) = k_b * x(3); % D + P <- D:P
lambda(3) = k_R * x(1); % D -> D + MR
lambda(4) = k_A * x(3); % D:P -> D:P + MA
lambda(5) = k_Rd * x(4); % MR -> 0
lambda(6) = k_Ad * x(5); % MA -> 0
end