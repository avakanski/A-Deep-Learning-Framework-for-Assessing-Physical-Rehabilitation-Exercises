% Codes for reading the data for the Deep Squat exercise
% from UI-PRMD dataset and preparing the data for NN processing
%
% The dataset can be downloaded from: https://www.webpages.uidaho.edu/ui-prmd/
% The data used in the paper can be found in the "Reduced Data Set" section
% 
%
% The outputs are the csv files Data_Correct and Data_Incorrect

clear; close;clc;

%% Load the files (Part 1)

% Number of subjects
nSubjects = 10;

% Number of repetitions per subject
nEpisodes = 10;

% Load the correct angle movements for DeepSquat
% To read text files use: M = csvread('m01_s01_e02_angles.txt');
for i=1:nSubjects-1
    for j=1:nEpisodes-1
        eval(['Train1_ang_cor{',int2str(i),',',int2str(j),'} = csvread(''Segmented Movements\Vicon\Angles\m01_s0',int2str(i),'_e0',int2str(j),'_angles.txt'');']);
    end
end

for i=1:nSubjects-1
    for j=nEpisodes
        eval(['Train1_ang_cor{',int2str(i),',',int2str(j),'} = csvread(''Segmented Movements\Vicon\Angles/m01_s0',int2str(i),'_e',int2str(j),'_angles.txt'');']);
    end
end

for i=nSubjects
    for j=1:nEpisodes-1
        eval(['Train1_ang_cor{',int2str(i),',',int2str(j),'} = csvread(''Segmented Movements\Vicon\Angles\m01_s',int2str(i),'_e0',int2str(j),'_angles.txt'');']);
    end
end

for i=nSubjects
    for j=nEpisodes
        eval(['Train1_ang_cor{',int2str(i),',',int2str(j),'} = csvread(''Segmented Movements\Vicon\Angles/m01_s',int2str(i),'_e',int2str(j),'_angles.txt'');']);
    end
end

% Load the incorrect angle movements for DeepSquat
% To read text files use: M = csvread('m01_s01_e02_angles_inc.txt');
for i=1:nSubjects-1
    for j=1:nEpisodes-1
        eval(['Test1_ang_inc{',int2str(i),',',int2str(j),'} = csvread(''Incorrect Segmented Movements\Vicon\Angles/m01_s0',int2str(i),'_e0',int2str(j),'_angles_inc.txt'');']);
    end
end

for i=1:nSubjects-1
    for j=nEpisodes
        eval(['Test1_ang_inc{',int2str(i),',',int2str(j),'} = csvread(''Incorrect Segmented Movements\Vicon\Angles\m01_s0',int2str(i),'_e',int2str(j),'_angles_inc.txt'');']);
    end
end

for i=nSubjects
    for j=1:nEpisodes-1
        eval(['Test1_ang_inc{',int2str(i),',',int2str(j),'} = csvread(''Incorrect Segmented Movements\Vicon\Angles/m01_s',int2str(i),'_e0',int2str(j),'_angles_inc.txt'');']);
    end
end

for i=nSubjects
    for j=nEpisodes
        eval(['Test1_ang_inc{',int2str(i),',',int2str(j),'} = csvread(''Incorrect Segmented Movements\Vicon\Angles/m01_s',int2str(i),'_e',int2str(j),'_angles_inc.txt'');']);
    end
end

%% Reorganize the data into an initial 1x100 array (Part 2)

% Correct sequences
k = 1;
for i=1:nSubjects
    for j=1:nEpisodes
        Correct_Ini{k} = Train1_ang_cor{i,j};
        k = k+1;
    end
end

% Total number of episodes
N = length(Correct_Ini);

% Dimensionality of sequences
nDim = size(Correct_Ini{1},2);

% Incorrect sequences
k = 1;
for i=1:nSubjects
    for j=1:nEpisodes
        Incorrect_Ini{k} = Test1_ang_inc{i,j};
        k = k+1;
    end
end

% Remove the first episodes for most of the subjects due to missing
% values during the data recording
% and also remove a few other inconsistent episodes
Correct_Red = Correct_Ini([2:10 12:20 22:30 32:40 42:50 51:60 62:63 65:70 72:80 82:84 86:90 91:100]);
Incorrect_Red = Incorrect_Ini([2:10 12:20 22:30 32:40 42:50 51:60 62:63 65:70 72:80 82:84 86:90 91:100]);

%% Perform linear alignment (Part 3)

% Number of correct sequences
nSequences = length(Correct_Red);

% Find the median sequence in length
for i = 1:nSequences 
    Len_Tr(i) = size(Correct_Red{i},1);
end
% Length_mean = ceil(mean(Len_Tr));
Length_mean = 240;

% Linear alignment to change the length of correct sequences to Length_mean
for i=1:nSequences 
    for j=1:nDim
        Correct_Aligned{i}(:,j)=interp1([1:size(Correct_Red{i},1)],Correct_Red{i}(:,j),linspace(1,size(Correct_Red{i},1),Length_mean));
    end
end

% Incorrect sequences
for i=1:nSequences 
    for j=1:nDim
        Incorrect_Aligned{i}(:,j)=interp1([1:size(Incorrect_Red{i},1)],Incorrect_Red{i}(:,j),linspace(1,size(Incorrect_Red{i},1),Length_mean));
    end
end

%% Transform into matrices with nSequences x nDim rows and sequence length columns(Part 4)

Correct_Xm = [Correct_Aligned{1}'];
for i = 2:nSequences 
    Correct_Xm = [Correct_Xm; Correct_Aligned{i}'];
end

Incorrect_Xm = [Incorrect_Aligned{1}'];
for i = 2:nSequences 
    Incorrect_Xm = [Incorrect_Xm; Incorrect_Aligned{i}'];
end

%% Re-center the data to obtain zero mean (Part 5)

Data_mean = repmat(mean(Correct_Xm,2), 1, size(Correct_Xm,2));
centered_Correct_Data = Correct_Xm - Data_mean;

% Recenter the incorrect sequences
Data_mean_inc = repmat(mean(Incorrect_Xm,2), 1, size(Incorrect_Xm,2));
centered_Incorrect_Data = Incorrect_Xm - Data_mean_inc;

% Scale the data between -1 and 1
scaling_value = ceil(max(max(max(centered_Correct_Data)),abs(min(min(centered_Correct_Data)))));
Data_Correct = centered_Correct_Data/scaling_value;
Data_Incorrect = centered_Incorrect_Data/scaling_value;

%% Save the data (Part 6)
csvwrite('Data_Correct.csv',Data_Correct)
csvwrite('Data_Incorrect.csv',Data_Incorrect)

