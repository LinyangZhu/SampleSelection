clc;
clear;
close all;

file_path = './';
file_name = 'SourceData.dat';
fid0 = fopen([file_path,'info_SelectData.dat'], 'w+');

coeff_std = 0.9999296;  
fprintf(fid0,'The critical cc is: %e \n', coeff_std);
m = importdata([file_path,file_name]);
row_sample = length(m);
fprintf('The number of original data: %d\n', row_sample);
fprintf(fid0,'The number of original data: %d \n',row_sample);

cri_var = 5;  % The critical value used to remove the outliers 
no_delete_total = 0;
for k=1:1
    no_delete = 0;
    data_norm = zscore(m(:,4:end));
    for i = 1:length(data_norm)
        for j = 1:size(data_norm,2)  
            if data_norm(i,j)>cri_var
                m(i-no_delete,:)=[];
                no_delete = no_delete+1;
                break
            end
        end
    end
    no_delete_total = no_delete_total + no_delete;
end
fprintf(fid0,'根据标准差删除的样本数为：%d \n',no_delete_total);
% Linearize the inputs and output
[row_m, col_m] = size(m);
input_output = m(1:row_m,4:col_m);  
[row_sample, col_sample] = size(input_output);
norm_input_output = linear_normalization(input_output);
%  Start the loop of Recursive Method
start_time=cputime;  
sample_NO = 1;
sample(1,:) = norm_input_output(1,:);
sample_selected(1) = 1;
for i=2:row_sample
    flag_1 = 0;
    coeff_max = 0;
    for k=1:sample_NO
        coeff = dot(norm_input_output(i,:),sample(k,:)) ...
            /(norm(norm_input_output(i,:))*norm(sample(k,:)));
        if coeff>coeff_std
            flag_1 = 1;
            break; 
        else
            coeff_max = max(coeff_max,coeff);
        end
    end
    if flag_1 == 1; continue; end
    if coeff_max < coeff_std
        sample_NO = sample_NO+1;
        sample(sample_NO,:) = norm_input_output(i,:);
%       sample_selected is the cell number of  sample_NOth sample in the ith sub-region
        sample_selected(sample_NO) = i;  
    end
end
time_cost = cputime - start_time;
fprintf('Time consumption: %fs \n',time_cost)
fprintf(fid0,'Time consumption: %fs \n',time_cost);
fprintf('The sample number afer sample selection is: %d\n', sample_NO);
fprintf(fid0,'The sample number afer sample selection is: %d \n',sample_NO);

% Save the datapairs after sample selection.(X_coor,Y_coor,q1,q2,q3,q4,output) 
file_save_name = 'data_select_Recursion.dat';
fid = fopen([file_path,file_save_name], 'w');
for i=1:sample_NO 
    for k=2:col_m
        fprintf(fid,'%e   ', m(sample_selected(i),k));
    end
    fprintf(fid,'\r\n');
end
fclose(fid);
status = fclose('all');

function out = linear_normalization(in)
[row,col] = size(in);
out = zeros(row,col);
for i=1:col
    max_v = max(in(:,i));
    min_v = min(in(:,i));
    out(:,i) = (in(:,i)-min_v)/(max_v-min_v);
end
end