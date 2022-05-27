clc;
clear;
close all;

file_path = './';
file_name = 'SourceData.dat'; % the data form is (yplus,X_coor,Y_coor,q1,q2,q3,q4,output)
fid0 = fopen([file_path,'Info_SelectData.dat'], 'w+');

coeff_std = 0.9999075;  
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
fprintf(fid0,'The NO of outliers is: %d \n',no_delete_total);
% Linearize the inputs and output
[row_m, col_m] = size(m);
input_output = m(1:row_m,4:col_m); 
[row_sample, col_sample] = size(input_output);
norm_input_output = linear_normalization(input_output);
%  Divide the flowfield into some sub-regions 
yplus = m(:,1);
lnyplus = log(yplus); % lny+
delta = 0.5; % Set the numerical interval
[max_lnyplus,~]=max(lnyplus);
cluster_NO = ceil(max_lnyplus/delta);  % The NO of sub-regions
CellNO_cluster=zeros(cluster_NO,1);  % The NO of datapairs in the sub-regions 
MaxCellNO_cluster = 1e3;
% The first column of cluster is the sequentional NO of the sub-region. 
% The second column of cluster is the cell NO 
cluster=zeros(cluster_NO,MaxCellNO_cluster); 
for i=1:row_m
    if lnyplus(i)<0; continue; end
    cluster_item = floor((lnyplus(i))/delta)+1;
    CellNO_cluster(cluster_item)=CellNO_cluster(cluster_item)+1;
    cluster(cluster_item,CellNO_cluster(cluster_item))= i; %i ±àºÅ
end
% Remove the sub-regions with zero datapairs
no_delete = 0;
for i=1:cluster_NO
    if cluster(i-no_delete,1) == 0
        cluster(i-no_delete,:) = [];
        CellNO_cluster(i-no_delete,:) = [];
        no_delete = no_delete +1;
    end
end
[row_NO, col_NO] = size(cluster);
fprintf('The number of sub-region is: %d \n', row_NO);
fprintf(fid0,'The number of sub-region is: %d \n', row_NO);
%  Start the loop of Physics-assisted Recursive Method
start_time=cputime; 
sample_NO_max = 1;
sample_NO_min = 1e3;
sample_selected_NO = 0;
sample_NO_cluster=zeros(cluster_NO,1);  % The datapair numbers in the sub-regions 
for i=1:row_NO
    flag_2 = 0;
    sample_NO = 1;
    sample(1,:) = norm_input_output(cluster(i,1),:);
    sample_selected(i,1) = cluster(i,1);
    for j=sample_NO+1:CellNO_cluster(i)
        if cluster(i,j)==0
            flag_2 = 1;
            break; 
        end
        flag_1 = 0;
        coeff_max = 0;
        for k=1:sample_NO
            coeff = dot(norm_input_output(cluster(i,j),:),sample(k,:)) ...
                /(norm(norm_input_output(cluster(i,j),:))*norm(sample(k,:)));
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
            sample(sample_NO,:) = norm_input_output(cluster(i,j),:);
%           sample_selected is the cell number of  sample_NOth sample in the ith sub-region
            sample_selected(i,sample_NO) = cluster(i,j); 
        end
    end
    if flag_2 == 1; continue; end
    sample_NO_max = max(sample_NO_max, sample_NO);
    sample_NO_min = min(sample_NO_min, sample_NO);
    sample_selected_NO = sample_selected_NO+sample_NO;
    sample_NO_cluster(i) =  sample_NO;
end
time_cost = cputime - start_time;
fprintf('Time consumption: %fs \n',time_cost)
fprintf(fid0,'Time consumption: %fs \n',time_cost);
fprintf('The sample number afer sample selection is: %d\n', sample_selected_NO);
fprintf(fid0,'The sample number afer sample selection is: %d \n',sample_selected_NO);
fclose(fid0);
cluster_sequence = 1:cluster_NO;
data_NO_comp = [transpose(cluster_sequence),CellNO_cluster,sample_NO_cluster];
% Save the datapair numbers before and after sample selection
fid2 = fopen([file_path,'data_NO_comp.dat'], 'w+');
fprintf(fid2,'%s \r\n','variables ="sub-region_sequence","datapair_NO_before_SS","datapair_NO_after_SS"');
for i=1:row_NO 
    fprintf(fid2,'%d   ',data_NO_comp(i,1));
    fprintf(fid2,'%d  %d \r\n', data_NO_comp(i,2),data_NO_comp(i,3));
end
fclose(fid2);

% Save the datapairs before sample selection according to subregion sequence
origin_data_no = 0;
file_save_name = 'data_origin.dat';
fid1 = fopen([file_path,file_save_name], 'w');
fprintf(fid1,'%s','title = "The original data"');
fprintf(fid1,'\r\n');
fprintf(fid1,'%s','variables ="NO", "q1","q2","q3","q4","output"');
fprintf(fid1,'\r\n');
for i=1:row_NO  
    yplus_lower = exp((i-1)*delta+1); yplus_upper = exp((i)*delta+1);
    fprintf(fid1,['zone T = ','"yplus=',num2str(yplus_lower),'-',num2str(yplus_upper),'"']);
    fprintf(fid1,'\r\n');
    for j=1:CellNO_cluster(i)
        fprintf(fid1,'%d   ',cluster(i,j));
        for k=1:col_sample
            fprintf(fid1,'%f   ', input_output(cluster(i,j),k));
        end
        origin_data_no = origin_data_no+1;
        fprintf(fid1,'\r\n');
    end
end
fclose(fid1);
% Save the datapairs after sample selection according to subregion sequence
selected_data_no = 0;
file_save_name = 'data_selected.dat';
fid = fopen([file_path,file_save_name], 'w');
fprintf(fid,'%s','title = "The selected samples"');
fprintf(fid,'\r\n');
fprintf(fid,'%s','variables ="NO", "q1","q2","q3","q4","output"');
fprintf(fid,'\r\n');
for i=1:row_NO 
    yplus_lower = exp((i-1)*delta+1); yplus_upper = exp((i)*delta+1);
    fprintf(fid,['zone T = ','"yplus=',num2str(yplus_lower),'-',num2str(yplus_upper),'"']);
    fprintf(fid,'\r\n');
    for j=1:sample_NO_max
        if sample_selected(i,j)==0;break;end
        fprintf(fid,'%d   ',sample_selected(i,j));
        for k=1:col_sample
            fprintf(fid,'%f   ', input_output(sample_selected(i,j),k));
        end
        selected_data_no = selected_data_no+1;
        fprintf(fid,'\r\n');
    end
end
fclose(fid);
% Save the datapairs after sample selection.(X_coor,Y_coor,q1,q2,q3,q4,output) 
file_save_name = 'data_selected_by_NewRecursion.dat';
fid = fopen([file_path,file_save_name], 'w');
for i=1:row_NO  
    for j=1:sample_NO_max
        if sample_selected(i,j)==0;break;end
        for k=2:col_m
            fprintf(fid,'%e   ', m(sample_selected(i,j),k));
        end
        fprintf(fid,'\r\n');
    end
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