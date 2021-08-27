function Output = Gillespie_Reduction(Stoi, propensity_functions, x_init, tspan, err_tol, rev_idx)
%% Argument Check
if nargin < 5
   error('Invalid Input')
end

if nargin == 6
   err_tol = 0.1; % Default error tolerance
end

%% Reduction
fprintf('Pre-processing...');
[var_num, react_num] = size(Stoi); % number of variables

stoi_tmp = Stoi;
Initial_temp = x_init;
RB_num = size(rev_idx, 1);
if RB_num > 0
    AT_thres = zeros(RB_num,1);
    L = zeros(RB_num,1);
    free = zeros(RB_num,2);
    cpx = zeros(RB_num,1);
    Kd = zeros(RB_num,1);
    
    for i = 1:RB_num
        idx1 = find(Stoi(:,rev_idx(i,1)) == 1);
        idx2 = find(Stoi(:,rev_idx(i,1)) == -1);
        if length(idx1) == 1
            free(i,:) = idx2;
            cpx(i) = idx1;
        else
            free(i,:) = idx1;
            cpx(i) = idx2;
        end
        input = zeros(var_num,1);
        input(free(i,1)) = 1;
        input(free(i,2)) = 1;
        input(cpx(i)) = 1;
        propensity = propensity_functions(input);
        
        Kd(i) = propensity(rev_idx(i,2))/propensity(rev_idx(i,1));  % Normalized Dissociation Const
        [AT_thres(i), L(i)] = QSSA_Threshold(Kd(i), err_tol);               % find threshord for the stQSSA and the slQSSA
        
        % Update Stoiciometric Matrix
        stoi_tmp(free(i,1),:) = Stoi(free(i,1),:) + Stoi(cpx(i),:);
        stoi_tmp(free(i,2),:) = Stoi(free(i,2),:) + Stoi(cpx(i),:);
        stoi_tmp(cpx(i,1),:) = zeros(1,react_num);
        
        % Update initial state
        Initial_temp(free(i,1)) = x_init(free(i,1)) + x_init(cpx(i));
        Initial_temp(free(i,2)) = x_init(free(i,2)) + x_init(cpx(i));
        Initial_temp(cpx(i,1)) = 0;
    end
end
fprintf('end\n');

%% Simulation
fprintf('Gillespie Simulation...');
rng('shuffle')
domain_len = length(tspan);
Tmax = tspan(end);
t = 0; %current time
k = Initial_temp; %current state
X = zeros(var_num, domain_len); % state record
j = 1;
%% iteration
while t <= Tmax
    % Calculate Propensity
    input = k;
    for i = 1:RB_num
        AT = input(free(i,1));
        BT = input(free(i,2));
        
        if  min(AT,BT) >= AT_thres
            input(free(i,1)) = (AT - BT - Kd(i) + sqrt((AT - BT - Kd(i))^2 + 4*AT*Kd))/2;
            input(cpx(i)) =  AT - input(free(i,1));
            input(free(i,1)) =  BT - input(cpx(i));
        else
            if AT <= BT
                input(free(i,1)) = LQSSA(AT,BT,Kd(i),L(i));
                input(cpx(i)) = AT - input(free(i,1));
                input(free(i,2)) =  BT - input(cpx(i));
            else
                input(free(i,2)) = LQSSA(BT,AT,Kd(i),L(i));
                input(cpx(i)) = BT - input(free(i,1));
                input(free(i,1)) =  AT - input(cpx(i));
            end
        end
    end
    
    propensity = propensity_functions(input);
    
    for i = 1:RB_num
        propensity(rev_idx(i,:)) = 0;
    end
    
    lambda = sum(propensity);
    if lambda == 0 %reactions do not happen
        while j <= scale_size
            X(:,j) = k;
            j = j + 1;
        end
        break
    end
    
    r = rand([2 1]); %two random numbers r1, r2
    T = -1/lambda*log(r(1));
    
    if t + T > Tmax %end condition
        while j <= domain_len
            X(:,j) = k;
            j = j + 1;
        end
        break
    end
    
    %choose the reaction
    propensity_sum = cumsum(propensity);
    reaction_index = find(propensity_sum > r(2) * lambda,1);
    
    %record the X(t)
    while tspan(j) < t + T
        X(:,j) = k;
        j = j + 1;
    end
    
    %update k and t
    k = k + stoi_tmp(:, reaction_index);
    t = t + T;
end
Output = X;
fprintf('end\n');
end