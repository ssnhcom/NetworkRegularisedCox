%% Divide sample into 2 groups: high risk and low risk
function [HighRisk, LowRisk] = getSampleGroups(data, time, status, beta)
    PI  = data' * beta;
    [orderPI, PIindex] = sort(PI,'descend');    
    M = length(PI); % number of patients
    N = int16(M * 0.4); % 40per of patients
    hindex = PIindex(1:N); % high risk index
    lindex  =  PIindex(M-N:M); % low risk index
    for i = 1:N
        HighRisk(i,1) = time(hindex(i));
        HighRisk(i,2) = status(hindex(i));
        
        LowRisk(i,1) = time(lindex(i));
        LowRisk(i,2) = status(lindex(i));
    end
end