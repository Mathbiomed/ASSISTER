% Developer: Yun Min Song
% Input: 
%   - Kd: Normalized dissociation constant (dissociation constant * volume)
%   - epsilon: Error tolerance 
% Output:
%   - AT_Bound: The minimum state of AT where the stQSSA approximate the
%   stochastic QSSA within a relative error of epsilon for given Kd value.
%   (i.e. when AT >= AT_thres, the relative error of the stQSSA to the 
%   stochastic QSSA is less than epsilon)
%   - L: The minimum number of states that the slQSSA requires to 
%   approximate the stochastic QSSA of A (<A>) within epsilon relative 
%   error when AT < AT_thres for given Kd value.
%   (i.e. when AT < AT_thres, the relative error of the L-state slQSSA to the 
%   stochastic QSSA is less than epsilon)

function [AT_thres, L] = QSSA_Threshold(Kd, epsilon)
%% Binary Search
L = 0;
R = ceil(1/(epsilon^2 * Kd)); % Upper bound of the stQSSA relative error
while L < R - 1
    m = floor((L+R) / 2);
    if tQerror(Kd,m) < epsilon
        R = m;
    else
        L = m + 1;
    end
end

AT_thres = R;
%% find minimum state number
A = LQSSA(AT_thres-1,AT_thres-1,Kd,AT_thres);
RA_L = Inf;
L = 1;
while RA_L > epsilon
    L = L+1;
    A_lq = LQSSA(AT_thres-1,AT_thres-1,Kd,L);
    RA_L = abs(A_lq/A - 1);
end

end

function RA_max = tQerror(Kd, AT)
RA_max = 0;
BT = AT;
while 1
    A = LQSSA(AT,BT,Kd,AT+1);
    A_tq = (AT-BT-Kd + sqrt((AT-BT-Kd)^2 + 4*AT*Kd))/2;
    if abs(A_tq/A - 1) >= RA_max
        RA_max  = abs(A_tq/A - 1);
    else
        break
    end
    BT = BT + 1;
end

end