% DTW_wthn_subjects: calculates DTW on the raw data for the within-subject case

clear; clc; close;

%% Load the data

% Correct repetitions
Data_NN = csvread('../../Data/Data_Correct.csv');

% Incorrect repetitions
Data_NN_inc = csvread('../../Data/Data_Incorrect.csv');

% Data dimensions
n_dim = 117;
% Number of timesteps
T1 = size(Data_NN,2);
% Number of repetitions
T2 = size(Data_NN,1)/n_dim;

% Tranform the data into cells
% Correct repetitions
Train_Data = cell(1,90);
for i=1:T2
    Train_Data{1,i} = Data_NN((i-1)*n_dim+1:i*n_dim,:)';
end

% Incorrect repetitions
Test_Data = cell(1,90);
for i=1:T2
    Test_Data{1,i} = Data_NN_inc((i-1)*n_dim+1:i*n_dim,:)';
end

% repetions for each subject
rt=[9 9 9 9 9 10 8 9 8 10];

%% Calculate the DTW distances

% Correct sequences
dtw_train = zeros(1,T2);
M=0;
for s=1:10
    for i=1:rt(1,s)
        for j=1:rt(1,s)
            dtw_train(M+i)=dtw_train(M+i)+dtw(Train_Data{M+i},...
                               Train_Data{M+j});
        end
        % scale data
        if rt(1,s)>0
           dtw_train(M+i)=dtw_train(M+i)/rt(1,s);
        else
           dtw_train(M+i)=dtw_train(M+i);
        end
    end
    M=M+rt(1,s);
end

% Incorrect sequences
dtw_test = zeros(1,T2);
N=0;
for s=1:10
    for i=1:rt(1,s)
        for j=1:rt(1,s)
            dtw_test(N+i) = dtw_test(N+i)+dtw(Test_Data{N+i},...
                               Train_Data{N+j});
        end
    % scale data
        if rt(1,s)>0
           dtw_test(N+i)=dtw_test(N+i)/rt(1,s);
        else
           dtw_test(N+i)=dtw_test(N+i);
        end
    end
    N=N+rt(1,s);
end

%% Scale and plot data

% Scale data
dtw_train = dtw_train/T1;
dtw_test = dtw_test/T1;

% Scale data in the [1,20] range
MAX = max(max(dtw_train),max(dtw_test));
MIN = min(min(dtw_train),min(dtw_test));
for s = 1: T2
    dtw_train(s) = 19*(dtw_train(s)-MIN)/(MAX-MIN)+1;
    dtw_test(s) = 19*(dtw_test(s)-MIN)/(MAX-MIN)+1;
end

% Plot
h = figure; 
plot(dtw_train,'go','LineWidth',2); hold on, ...
plot(dtw_test,'rs','LineWidth',2);
xlabel('Sequence Number', 'fontsize',18);
ylabel('DTW Distance', 'fontsize',18);
% Title: DTW Distance Within Subjects - No Dimensionality Reduction
title('DTW Distance WS NDR','fontsize',18);
legend({'Correct Sequences','Incorrect Sequences'}, 'fontsize',16,'location','NW')
set(gca,'box','off','fontweight','bold','LineWidth',2);
set(gcf,'Units','inches','position',[0 0 5.5 4.5]);
% print(h,'../../Results/DTW_Distance_WS_NDR','-dpng','-r300');

%% Separation degree betwween the correct and incorrect sequences
SD = zeros(1,1);
for i=1:T2
    for j=1:T2
        SD=SD+(dtw_test(i)-dtw_train(j))/(abs(dtw_test(i))+abs(dtw_train(j)));
    end
end
SD=SD/T2/T2;

