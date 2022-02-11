function [solutions] = linearspline(t1, v1, t2)
% Linear spline to interpolate the values. 
% No extrapolation is allowed
% Given data: (t1,v1), return data interpolated at t2
% assume evenly spaced intervals in both t1 and t2.

n1 = length(t1);
n2 = length(t2);

dt1 = t1(2)-t1(1);

n = n2;

count = 0;

if(t2(end)<=t1(end))
    n = n2;
    t = t2;
else
    for i=1:n2
        t2Val = t2(i);
        if(t2Val <= t1(end))
            count = count+1;
            t(count) = t2Val;
        end
    end
    n = count;
end

v = zeros(1,n);
diffN = (v1(2:end)-v1(1:end-1))./dt1;

for i=1:n
    t2Val = t(i);
    tempVal = 0;
   for j=2:n1
         t1p1 = t1(j);
         t1Val = t1(j-1);
         if(t2Val >= t1Val && t2Val <= t1p1)
             tempVal = v1(j-1)+diffN(j-1)*(t2Val-t1Val);
             break;
         end
   end
   if tempVal < 0   %ask about this
           v(i) = 0;
   else
           v(i) = tempVal;
   end
end

solutions = zeros(1,15);

for i = 1:length(v)
    
    solutions(i) = v(i);
end

