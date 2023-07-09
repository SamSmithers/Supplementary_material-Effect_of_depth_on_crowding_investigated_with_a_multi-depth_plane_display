% This is a copy of the MATLAB script used to accumulate the raw data from each experiment 
% and calculate the perceptual error as reported in Smithers et al.- "Large differences in 
% target-flanker depth can increase crowding- evidence from a multi-depth plane display".  

% Script written by Dr Samuel P. Smithers, Northeastern University, 2022-2023
% Last edited July 2023

% Corresponding authors SPS (s.smithers@northeastern.edu) and PJB (p.bex@northeastern.edu)

% The script first calculates the report error (difference between the reported gap 
% orientation and the actual gap orientation) for each individual stimulus presentation. 
% It then calculates the circular standard deviation for each stimulus condition (based 
% on 10 repeats (8 in Exp 4)) and saves this output as an accumulated results file. This 
% script also has the option to do this using only the repeats in which the subject 
% reported seeing the target inside the flanker ring (as reported in the supplementary 
% material).

% This script uses the circular statistics toolbox:
% https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics

clear
close all
SaveData = "yes"; % If "yes" the script will save the results output as a .csv file. 
IncludeTargetInsideRingOnly = "no";
% ^^^ If IncludeTargetInsideRingOnly = "no": circular SD for each stimulus condition is calculated from all 10 repeats (8 in Exp 4).
% ^^^ If IncludeTargetInsideRingOnly = "yes": circular SD for each stimulus condition is calculated based only on trials in which
% the observer reported seeing the target inside the flanker ring.

%% Choose which experiment to process by uncommenting the appropriate 'Savename' and 'DirectoryPath' below.
% Exp 1
SaveName = 'Exp1_accumulatedPerErr';
DirectoryPath = 'PATH WAY TO FOLDER CONTAINING "Exp 1 raw data"\'; % Folder must contain only the raw data from Experiment 1.

% Exp 2
% SaveName = 'Exp2_accumulatedPerErr';
% DirectoryPath = 'PATH WAY TO FOLDER CONTAINING "Exp 2 raw data"\'; % Folder must contain only the raw data from Experiment 2.

% Exp 3
% SaveName = 'Exp3_accumulatedPerErr';
% DirectoryPath = 'PATH WAY TO FOLDER CONTAINING "Exp 3 raw data"\'; % Folder must contain only the raw data from Experiment 3.

% Exp 4
% SaveName = 'Exp4_accumulatedPerErr';
% DirectoryPath = 'PATH WAY TO FOLDER CONTAINING "Exp 4 raw data"\'; % Folder must contain only the raw data from Experiment 4.

% Exp 5
% SaveName = 'Exp5_accumulatedPerErr';
% DirectoryPath = 'PATH WAY TO FOLDER CONTAINING "Exp 5 raw data"\'; % Folder must contain only the raw data from Experiment 5.

% DATA FROM EXPERIENCED SUBJECTS
% Exp 1
% SaveName = 'Exp1_accumulatedPerErr_experiencedSubjects';
% DirectoryPath = 'PATH WAY TO FOLDER CONTAINING "Exp 1 raw data (experienced subjects)"\'; % Folder must contain only the raw data from experienced subjects for Experiment 1.

% Exp 2
% SaveName = 'Exp2_accumulatedPerErr_experiencedSubjects';
% DirectoryPath = 'PATH WAY TO FOLDER CONTAINING "Exp 2 raw data (experienced subjects)"\'; % Folder must contain only the raw data from experienced subjects for Experiment 2.

% Exp 3
% SaveName = 'Exp3_accumulatedPerErr_experiencedSubjects';
% DirectoryPath = 'PATH WAY TO FOLDER CONTAINING "Exp 3 raw data (experienced subjects)"\'; % Folder must contain only the raw data from experienced subjects for Experiment 3.

% Exp 4
% SaveName = 'Exp4_accumulatedPerErr_experiencedSubjects';
% DirectoryPath = 'PATH WAY TO FOLDER CONTAINING "Exp 4 raw data (experienced subjects)"\'; % Folder must contain only the raw data from experienced subjects for Experiment 4.

% Exp 5
% SaveName = 'Exp5_accumulatedPerErr_experiencedSubjects';
% DirectoryPath = 'PATH WAY TO FOLDER CONTAINING "Exp 5 raw data (experienced subjects)"\'; % Folder must contain only the raw data from experienced subjects for Experiment 5.


%% Calculate perceptual error (i.e. circular standard deviation) for each stimulus condition and accumulate results into a single file. 
Files = dir(fullfile(DirectoryPath, '*_rawData.csv')); %Get list of all raw data files within the folder.

CollatedProcessedData = table();

% Loop through all files in the folder.
for FileNo = 1:length(Files)
    FileName = Files(FileNo).name;
    fullFileName = fullfile(DirectoryPath, FileName);
    disp('Now processing')
    disp(FileName)

    % Get data from file.
    data = readtable(fullFileName);

    % Calculate report error for each row.
    for ii = 1:height(data)
        Actual_Ori_radian = deg2rad(data.Presented_Landolt_Ori_in_deg(ii));
        Reported_Ori_radian = deg2rad(data.Reported_Landolt_Ori_in_deg(ii));
        ang = angdiff(Actual_Ori_radian,Reported_Ori_radian);
        data.ReportError_radians(ii)=ang;
    end

    % Loop through the different stimulus conditions.
    ProcessedData = table(); % Make an empty table to fill.

    for ii = 1:max(data.Condition_number)
        subdata = data(data.Condition_number==ii,:);

        % Extract data and info needed for analysis.
        ProcessedData.ID(ii) = subdata.ID(1);
        ProcessedData.Condition_number(ii) = subdata.Condition_number(1);
        ProcessedData.Target_Screen(ii) = subdata.Target_Screen(1);
        ProcessedData.Flanker_Screen(ii) = subdata.Flanker_Screen(1);
        ProcessedData.Fixation_Screen(ii) = subdata.Fixation_Screen(1);
        ProcessedData.Target_Flanker_Spacing_in_deg(ii) = subdata.Target_Flanker_Spacing_in_deg(1);

        if IncludeTargetInsideRingOnly == "no"
            try
                % Count the number of times each perceived target location option was chosen.
                countdata = subdata(subdata.Perceived_TargetRing_Position== "target in center of ring",:);
                TPcount1 = height(countdata.Perceived_TargetRing_Position);
                ProcessedData.PerceivedPosCount_TCenter(ii) = TPcount1;
                countdata = subdata(subdata.Perceived_TargetRing_Position== "target not in center of ring",:);
                TPcount2 = height(countdata.Perceived_TargetRing_Position);
                ProcessedData.PerceivedPosCount_TOffCenter(ii) = TPcount2;
                countdata = subdata(subdata.Perceived_TargetRing_Position== "ring obstructs target",:);
                TPcount3 = height(countdata.Perceived_TargetRing_Position);
                ProcessedData.PerceivedPosCount_TObstructed(ii) = TPcount3;
                countdata = subdata(subdata.Perceived_TargetRing_Position== "target outside ring",:);
                TPcount4 = height(countdata.Perceived_TargetRing_Position);
                ProcessedData.PerceivedPosCount_TOutside(ii) = TPcount4;
                countdata = subdata(subdata.Perceived_TargetRing_Position== "unsure or no ring",:);
                TPcount5 = height(countdata.Perceived_TargetRing_Position);
                ProcessedData.PerceivedPosCount_noRingOrUnsure(ii) = TPcount5;
                TotalTPcount = TPcount1 + TPcount2 + TPcount3 + TPcount4 + TPcount5;
                % Check that the count total is correct.
                if TotalTPcount ~= 10 && TotalTPcount ~= 8
                    disp("*************************************************************")
                    disp("WARNING: The number of repeats is not 10 (or 8 in the case of Exp 4).")
                    disp("*************************************************************")
                    breakpoint() % The above warning indicates an issues with the data that should be investigated before continuing.
                end
            catch % We use a try and catch here because this step is not applicable to Experiment 5. 
            end
        elseif IncludeTargetInsideRingOnly == "yes" 
            % Exclude trials in which the subject did not perceive the target as being inside the flanker ring.
            try
                subdata = subdata( ...
                    subdata.Perceived_TargetRing_Position=="target in center of ring" | ...
                    subdata.Perceived_TargetRing_Position=="target not in center of ring" | ...
                    subdata.Flanker_Screen=="no flankers",:); % All data from the control condition is included.
            catch
                disp("*******************************************************************************************************************")
                disp("'Perceived_TargetRing_Position' is not applicable to Experiment 5. 'IncludeTargetInsideRingOnly' should be changed to 'no'.")
                disp("*******************************************************************************************************************")
                breakpoint()
            end
        end

        DataCount = height(subdata);
        
        if DataCount >2 % When IncludeTargetInsideRingOnly == "yes", circular SD is only calculated if there 
                        % are 3 or more data points for that condition.

            % Calculate the circular standard deviation.
            % circ_std() Computes circular standard deviation for circular data.
            %   Output:
            %     s     angular deviation
            %     C_SD  circular standard deviation
            [s, C_SD] = circ_std(subdata.ReportError_radians, [], [], []);
            ProcessedData.CircSD_deg(ii) = rad2deg(C_SD);

            if contains(SaveName,"_experiencedSubjects")
                % Calculate 95% bootstrapped confidence intervals for circular standard deviation
                BootStrappedSDs=zeros(1,1000);
                for bootNum = 1:1000 
                    SampleData = randsample(subdata.ReportError_radians,length(subdata.ReportError_radians),'true');
                    [s, C_SD] = circ_std(SampleData, [], [], []);
                    BootStrappedSDs(bootNum) = rad2deg(C_SD);
                end
                ProcessedData.BootStrappedCI(ii) = 1.96*std(BootStrappedSDs);
            end

        elseif DataCount <3
            ProcessedData.CircSD_deg(ii) = nan;
        end

    end

    % Collate all processed data.
    CollatedProcessedData = [CollatedProcessedData;ProcessedData];
end

% Save data file.
if SaveData == "yes"
    if IncludeTargetInsideRingOnly == "yes"
        name = sprintf('%s_TargetInsideRing.csv',SaveName);
    else
        name = sprintf('%s.csv',SaveName);
    end
    writetable(CollatedProcessedData,name,'Delimiter',',')
end
