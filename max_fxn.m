function [solutions] = max_fxn(str1,str2, bleb_start, frames, plotn, bleb_frame1, bleb_frame2)
% function to plot (if plotn = 1) the max values from cell data over time and return
% vector with the max values

% initialize t0
t0 = bleb_start;

% Preallocate space
solutions = zeros(1,15);

% LOAD CELL DATA 
bleb_frame = 1; % To extract part of your data begining with column 1
str = 'excel-sheets/'+str2;
bleb_MhcA_struct = importdata(str);
bleb_MhcA_data = bleb_MhcA_struct.data;
[temp1, ~] = myosin_cortex_fxn(frames, bleb_frame, bleb_MhcA_data, 0);
max_MhcA_vec_A = zeros(1,15);
bleb_frame = bleb_frame1; % To extract part of your data begining with column 18
[temp2, ~] = myosin_cortex_fxn(frames, bleb_frame, bleb_MhcA_data, 0);
max_MhcA_vec_B = zeros(1,15);
bleb_frame = bleb_frame2; % To extract part of your data begining with column 35
[temp3, T] = myosin_cortex_fxn(frames, bleb_frame, bleb_MhcA_data, 0);
max_MhcA_vec_C = zeros(1,15);


for i = t0:frames
     max_MhcA_vec_A(i-t0+1) = temp1(i);
     max_MhcA_vec_B(i-t0+1) = temp2(i);
     max_MhcA_vec_C(i-t0+1) = temp3(i);
end

% Find Average of 3 Values 
for i = 1:length(max_MhcA_vec_A)
    solutions(i) = (max_MhcA_vec_A(i) + max_MhcA_vec_B(i) + ... 
       max_MhcA_vec_C(i))./3;
end

%PLOT MAX VALUES OVER TIME
if plotn == 1  
plotMax(T(t0-1:frames), temp1, temp2, temp3, str1)
end

end

