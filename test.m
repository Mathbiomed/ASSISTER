Stoi = [-1 -1 1 0 0 ; 1 1 -1 0 0 ; 0 0 0 1 0;...
         0 0 0 0 1; 0 0 0 -1 0; 0 0 0 0 -1]'; %Stoichiometric

tspan = 0:200:1600;
err_tol = 0.1;
x_init = [10,10,0,0,0]';
rev_idx = [1 2];
prop_functions = @propensity_functions;

species = Gillespie_Reduction(gamma, prop_functions, x_init, tspan, err_tol, rev_idx);
disp(species)

function f = propensity_functions(k)
%fast reaction const
kf = 100; %binding
kb = 1; %unbinding

%slow reaction const
k1 = 0.1; %M1 transcription
k2 = 0.001; %M2 transcription
kd1 = 0.001; %M1 Degradation
kd2 = 0.001; %M2 Degradation

f = [kf*k(1)*k(2); kb*k(3); k1*k(1);
            k2*k(3); kd1*k(4); kd2*k(5)];
end