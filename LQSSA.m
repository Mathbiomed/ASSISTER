% Developer: Yun Min Song
% Input: 
%   - AT, BT: Total numbers of the molecules of the binding species.
%   - Kd: Normalized dissociation constant (dissociation constant x volume)
%   - N: Number of states used to calculate the slQSSA 
% Output: 
%   - f: The slQSSA of A in the reversible binding reaction A + B <=> C 
%   which is the approximation of of the stochastic QSSA of A (stationary 
%   average number of A molecules)

function f = LQSSA(AT,BT,Kd,N) 
    if (N < 1) || (AT < 0) || (BT < 0) || (Kd < 0)                         
        error('Error: Please check input parameters')
    end
    
    max_state = min(AT+1,N);                                                
    
    if max_state == 1
       f = 0; 
       return
    end
    numer = 0;
    denumer = 0;
    if AT < BT
        for l = (max_state-2):(-1):0
            g = Kd*(AT - l)/(l+1)/(BT - AT + l + 1);
            numer = (numer + l + 1)*g;
            denumer = (denumer + 1)*g;
        end
        f = numer/(denumer+1);
        return
    else
        for l = (max_state-2):(-1):0
            g = Kd*(BT - l)/(l+1)/(AT - BT + l + 1);
            numer = (numer + l + 1)*g;
            denumer = (denumer + 1)*g;
        end
        f = AT - BT + numer/(denumer+1);
        return
    end 
end