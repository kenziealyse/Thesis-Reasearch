function k5plus = kfiveplus(t, k5plusinit, a, tSteadyState)

%KFIVEPLUS Summary of this function goes here
%   Detailed explanation goes here

if (0 < t) && (t < tSteadyState)
    
    k5plus = k5plusinit;
    
else 
    
    k5plus = k5plusinit*exp(-a*t);
    
end

end

