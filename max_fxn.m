function [max_MhcA_vec_A, max_MhcA_vec_B, max_MhcA_vec_C] = max_fxn(str1,str2, bleb_start, frames, plotn, bleb_frame1, bleb_frame2)

t0 = bleb_start;

% LOAD CELL DATA 
bleb_frame = 1; % To extract part of your data begining with column 1
str = str1;
bleb_MhcA_struct = importdata(str2);
bleb_MhcA_data = bleb_MhcA_struct.data;


[temp1, ~] = myosin_cortex_fxn(frames, bleb_frame, bleb_MhcA_data, 0);

max_MhcA_vec_A = zeros(1,15);


% LOAD CELL DATA
bleb_frame = bleb_frame1; % To extract part of your data begining with column 18

[temp2, ~] = myosin_cortex_fxn(frames, bleb_frame, bleb_MhcA_data, 0);

max_MhcA_vec_B = zeros(1,15);


% LOAD CELL DATA
bleb_frame = bleb_frame2; % To extract part of your data begining with column 35

[temp3, T] = myosin_cortex_fxn(frames, bleb_frame, bleb_MhcA_data, 0);

max_MhcA_vec_C = zeros(1,15);

for i = t0:frames
    
    max_MhcA_vec_A(i) = temp1(i);
    
    max_MhcA_vec_B(i) = temp2(i);
    
    max_MhcA_vec_C(i) = temp3(i);
    
end


%PLOT MAX VALUES OVER TIME

if plotn == 1
    
plotMax(T(t0:end), temp1, temp2, temp3, str)

end

end

