% analyzing the Kinematics data and ECoG data. 
% in this version; save:
% 1- frequency of recording  
% 2- the index for pos max and vel max for later sepearting flexion and extension
% 3- the brain areas; layout for all good channels 
% 4- related processed ECoG signals

close all;
clear all;
clc;

%% General parameters that should be specified at the begining

% choosing brain area; could have effect on median referencing(?)
% options:
%Selected_Chs;  case 1 all electrodes
%Selected_PMChs; case 2 
%Selected_SMChs; case 3
%Selected_PMSMChs; case 4
%Selected_HandChs; case 5
Brain_part=1;

%type of referencing the ECoG
Reference=1; % for median
%Reference=2; % for mean 

%% bad channels and time stamps for trials
load('E:\ECoGLeapMotion\DataPatientTwo\ECoGData\ECoG_data.mat');

% Final list of bad channels by observation and other methods
BadChs=[65 66 128 128+1 128+2 128+16 128+18 128+20 128+24 128+30 ...
    128+31 128+32 128+44 128+62 128+64];

% make channels Ready for median or mean
All_Index=ones(size(lfp,1),1);
All_Index(BadChs')=0;

% only 256 Channels were used for real recording
Nch_Record=256;
Selected_Chs=logical(All_Index(1:Nch_Record,1));

% time markers in ecog
[Nch,N] = size(lfp);
ecog_time = (0:N-1)/Fs;
anin = anin(1,1:N);

% from Daniel function: start and end of movements
[trial_start_time,trial_end_time] = get_trial_times(anin,Fs);

trial_start=trial_start_time*Fs;
trial_stop=trial_end_time*Fs;

%% Finding the precentral (for Primary Motor=PM) and postcentral (somatosensory=SM) electrodes/channels

 load('E:\ECoGLeapMotion\DataPatientTwo\ImagingData\TDT_elecs_all.mat')
 
 PM_chs=[];
 SM_chs=[];
 PMSM_chs=[];
 
 for i=65:320
     
     % removing/passing bad channels 
     if ismember(i-64,BadChs)
         continue 
         
      % finding PM channels   
     elseif strcmp(anatomy{i,4},'precentral')
         PM_chs=[PM_chs;(i-64)];
         PMSM_chs=[PMSM_chs;(i-64)];
         
      % finding SM channels   
     elseif strcmp(anatomy{i,4},'postcentral')
         SM_chs=[SM_chs;(i-64)];
         PMSM_chs=[PMSM_chs;(i-64)];
         
     end   
     
     
 end
 

% logic PM channels
All_Index_PM=zeros(size(lfp,1),1);
All_Index_PM(PM_chs)=1;
Selected_PMChs=logical(All_Index_PM(1:Nch_Record,1));
 

% logic SM channels
All_Index_SM=zeros(size(lfp,1),1);
All_Index_SM(SM_chs)=1;
Selected_SMChs=logical(All_Index_SM(1:Nch_Record,1));

% logic PMSM channels
All_Index_PMSM=zeros(size(lfp,1),1);
All_Index_PMSM(PMSM_chs)=1;
Selected_PMSMChs=logical(All_Index_PMSM(1:Nch_Record,1));


%% Finding the hand area electrodes/channels with observation

% hand area channels
for i=0:4
    HandChs(((i*7)+1):((i*7)+7))=(122:128)+i*16;
end 

 All_Index_Hand=zeros(Nch_Record,1);
 for i=1:Nch_Record
     
     % removing/passing bad channels 
     if ismember(i,BadChs)
         continue 
         
      % finding hand channels   
     elseif ismember(i,HandChs)
         All_Index_Hand(i,1)=1;   
     end     
 end
 
 % logic Hand channels
Selected_HandChs=logical(All_Index_Hand);

%% Determining brain section for analysis e.g.: median ref

switch Brain_part
    
    case 1
        Brain_Area=Selected_Chs;
    case 2
        Brain_Area=Selected_PMChs;
    case 3
        Brain_Area=Selected_SMChs;
    case 4
        Brain_Area=Selected_PMSMChs;
    case 5
        Brain_Area=Selected_HandChs;
end

 %% All Kinematics data for all movements: loading/preparation 
num_trials =8;

for trial=1:num_trials
    % load trial data
    filename = 'E:\ECoGLeapMotion\DataPatientTwo/LeapMotionData/20180301/10h02m38s/EC171_20180301_10h02m38s_trial00';
    load([filename,num2str(trial),'.mat']);
   
    % grab kinematics
    [time_trial,hand_trial,fingers_trial] = leap_motion_kinematics(data.leap_frame);
   
    % time stamps for kinematics data
    time_AllTrials{trial} = time_trial; 
    % hand pos
    hand{trial}.palm_pos = hand_trial.palm_pos;
    hand{trial}.palm_vel = hand_trial.palm_vel;
    hand{trial}.palm_normal = hand_trial.palm_normal;
    % finger positions
    for f=1:5
        fingers{trial}(f).pos = fingers_trial(f).pos;
        fingers{trial}(f).vel = fingers_trial(f).vel;
        fingers{trial}(f).dir = fingers_trial(f).dir;
    end
end

% smooth all the trajectories and calc velocities
for trial=1:num_trials
    for f=1:5
        fingers{trial}(f).spos = smooth_kinematics(fingers{trial}(f).pos);
        fingers{trial}(f).svel = smooth_kinematics(fingers{trial}(f).vel);
        fingers{trial}(f).tvel = sqrt(sum(fingers{trial}(f).svel.^2,2)); % tangential
    end
end

%% 1- generating the max pos and max vel using the dominant dimension

%% 1-1 F1 movement analysis: Finding the peak points for transition from flextion to extention or extension to flextion

Xq_F1=0:1/(Fs):time_AllTrials{1,1}(end);
Posq_F1=interp1(time_AllTrials{1,1},fingers{1,1}(1).spos,Xq_F1);

half_window=200; %points; about 1 sec: for finding local max and local min in position data

% for finger 1 
[F1_n,F1_m]=size(Xq_F1');

% finding max and min time stamps
j=1;
jj=1;
for i=15*half_window:F1_n-(half_window+1)
    
    logic1=(Posq_F1(i-half_window:i+half_window,1)< Posq_F1(i,1));
    
    if sum(logic1)==2*half_window
        FingersKinInfo.Finger(1).Fs_MaxPos(j,1)=i;
        j=j+1;
    end
    
    logic2=(Posq_F1(i-half_window:i+half_window,1)> Posq_F1(i,1));
    
    if sum(logic2)==2*half_window
        FingersKinInfo.Finger(1).Es_MaxPos(jj,1)=i;
        jj=jj+1;
    end
    
    
end

% check on the plot
figure;
plot(Posq_F1(:,1))
hold on; plot(FingersKinInfo.Finger(1).Fs_MaxPos,-30,'og')
hold on; plot(FingersKinInfo.Finger(1).Es_MaxPos,-5,'or')

%% 1-2 F1 movement analysis: Finding the max velocity during transition from flextion to extention or extension to flextion

Velq_F1=interp1(time_AllTrials{1,1},fingers{1,1}(1).svel,Xq_F1);
AbsVel_F1=sqrt(Velq_F1(:,1).^2+Velq_F1(:,2).^2+Velq_F1(:,3).^2);

Index_F1=sort([1; FingersKinInfo.Finger(1).Fs_MaxPos;...
    FingersKinInfo.Finger(1).Es_MaxPos; F1_n]);


j=1;
for i=1:2:(length(Index_F1)-1)
    Section=AbsVel_F1(Index_F1(i):(Index_F1(i+1)));
    [m_S n_S]=max(Section);
    FingersKinInfo.Finger(1).Fs_MaxVel(1,j)=(n_S+Index_F1(i)-1);  
    j=j+1;
end

j=1;
for i=2:2:(length(Index_F1)-1)
    Section=AbsVel_F1(Index_F1(i):(Index_F1(i+1)));
    [m_S n_S]=max(Section);
    FingersKinInfo.Finger(1).Es_MaxVel(1,j)=(n_S+Index_F1(i)-1);
    j=j+1;
end


% check on the plot
figure;
plot(AbsVel_F1)
hold on; plot(FingersKinInfo.Finger(1).Fs_MaxVel,100,'ok')
hold on; plot(FingersKinInfo.Finger(1).Es_MaxVel,150,'ob')

figure;
plot(Posq_F1(:,1))
hold on; plot(FingersKinInfo.Finger(1).Fs_MaxPos,-30,'og')
hold on; plot(FingersKinInfo.Finger(1).Es_MaxPos,-5,'or')
hold on; vline(FingersKinInfo.Finger(1).Fs_MaxVel,'g--')
hold on; vline(FingersKinInfo.Finger(1).Es_MaxVel,'r--')

FingersKinData.Finger(1).PurePos=Posq_F1(:,1);
FingersKinData.Finger(1).PureVel=AbsVel_F1;


%% 1-3 F2 movement analysis: Finding the peak points for transition from flextion to extention or extension to flextion

Xq_F2=0:1/(Fs):time_AllTrials{1,2}(end);
Posq_F2=interp1(time_AllTrials{1,2},fingers{1,2}(2).spos,Xq_F2);

half_window=200; %points; about 1 sec: for finding local max and local min in position data

% for finger 1 
[F2_n,F2_m]=size(Xq_F2');



% finding max and min time stamps
j=1;
jj=1;

for i=10*half_window:F2_n-(half_window+1)
    
    logic1=(Posq_F2(i-half_window:i+half_window,3)< Posq_F2(i,3));
    
    if sum(logic1)==2*half_window
        FingersKinInfo.Finger(2).Fs_MaxPos(j,1)=i;
        j=j+1;   
    end
    
    logic2=(Posq_F2(i-half_window:i+half_window,3)> Posq_F2(i,3));
    
    if sum(logic2)==2*half_window
        FingersKinInfo.Finger(2).Es_MaxPos(jj,1)=i;
        jj=jj+1;  
    end
    
    
end

% check on the plot
figure;
plot(Posq_F2(:,3))
hold on; plot(FingersKinInfo.Finger(2).Fs_MaxPos,-80,'og')
hold on; plot(FingersKinInfo.Finger(2).Es_MaxPos,-50,'or')

%% 1-4 F2 movement analysis: Finding the max velocity during transition from flextion to extention or extension to flextion

Velq_F2=interp1(time_AllTrials{1,2},fingers{1,2}(2).svel,Xq_F2);
AbsVel_F2=sqrt(Velq_F2(:,1).^2+Velq_F2(:,2).^2+Velq_F2(:,3).^2);

Index_F2=sort([1; FingersKinInfo.Finger(2).Fs_MaxPos;...
    FingersKinInfo.Finger(2).Es_MaxPos; F2_n]);

j=1;
for i=1:2:(length(Index_F2)-1)
    Section=AbsVel_F2(Index_F2(i):(Index_F2(i+1)));
    [m_S n_S]=max(Section);
    FingersKinInfo.Finger(2).Fs_MaxVel(1,j)=(n_S+Index_F2(i)-1);
    j=j+1;
    
end

j=1;
for i=2:2:(length(Index_F2)-1)
    Section=AbsVel_F2(Index_F2(i):(Index_F2(i+1)));
    [m_S n_S]=max(Section);
    FingersKinInfo.Finger(2).Es_MaxVel(1,j)=(n_S+Index_F2(i)-1);
    j=j+1;
end


% check on the plot
figure;
plot(AbsVel_F2)
hold on; plot(FingersKinInfo.Finger(2).Fs_MaxVel,100,'ok')
hold on; plot(FingersKinInfo.Finger(2).Es_MaxVel,150,'ob')

figure;
plot(Posq_F2(:,3))
hold on; plot(FingersKinInfo.Finger(2).Fs_MaxPos,-80,'og')
hold on; plot(FingersKinInfo.Finger(2).Es_MaxPos,-50,'or')
hold on; vline(FingersKinInfo.Finger(2).Fs_MaxVel,'g--')
hold on; vline(FingersKinInfo.Finger(2).Es_MaxVel,'r--')

FingersKinData.Finger(2).PurePos=Posq_F2(:,3);
FingersKinData.Finger(2).PureVel=AbsVel_F2;

%% 1-5 F3 movement analysis: Finding the peak points for transition from flextion to extention or extension to flextion

Xq_F3=0:1/(Fs):time_AllTrials{1,3}(end);
Posq_F3=interp1(time_AllTrials{1,3},fingers{1,3}(3).spos,Xq_F3);

half_window=250; %points; about 1 sec: for finding local max and local min in position data

% for finger 1 
[F3_n,F3_m]=size(Xq_F3');

% finding max and min time stamps
j=1;
jj=1;
for i=10*half_window:F3_n-(half_window+1)
    
    logic1=(Posq_F3(i-half_window:i+half_window,3)< Posq_F3(i,3));
    
    if sum(logic1)==2*half_window
        FingersKinInfo.Finger(3).Fs_MaxPos(j,1)=i;
        j=j+1;  
    end
    
    logic2=(Posq_F3(i-half_window:i+half_window,3)> Posq_F3(i,3));
    
    if sum(logic2)==2*half_window
        FingersKinInfo.Finger(3).Es_MaxPos(jj,1)=i;
        jj=jj+1;
    end
    
    
end

FingersKinInfo.Finger(3).Fs_MaxPos=FingersKinInfo.Finger(3).Fs_MaxPos([1:3,5:end]);
% check on the plot
figure;
plot(Posq_F3(:,3))
hold on; plot(FingersKinInfo.Finger(3).Fs_MaxPos,-80,'og')
hold on; plot(FingersKinInfo.Finger(3).Es_MaxPos,-50,'or')

%% 1-6 F3 movement analysis: Finding the max velocity during transition from flextion to extention or extension to flextion

Velq_F3=interp1(time_AllTrials{1,3},fingers{1,3}(3).svel,Xq_F3);
AbsVel_F3=sqrt(Velq_F3(:,1).^2+Velq_F3(:,2).^2+Velq_F3(:,3).^2);

Index_F3=sort([1; FingersKinInfo.Finger(3).Fs_MaxPos;...
    FingersKinInfo.Finger(3).Es_MaxPos; F3_n]);

j=1;
for i=1:2:(length(Index_F3)-1)
    Section=AbsVel_F3(Index_F3(i):(Index_F3(i+1)));
    [m_S n_S]=max(Section);
    FingersKinInfo.Finger(3).Fs_MaxVel(1,j)=(n_S+Index_F3(i)-1);
    j=j+1;
end

j=1;
for i=2:2:(length(Index_F3)-1)
    Section=AbsVel_F3(Index_F3(i):(Index_F3(i+1)));
    [m_S n_S]=max(Section);
    FingersKinInfo.Finger(3).Es_MaxVel(1,j)=(n_S+Index_F3(i)-1);
    j=j+1;
end

% check on the plot
figure;
plot(AbsVel_F3)
hold on; plot(FingersKinInfo.Finger(3).Fs_MaxVel,100,'ok')
hold on; plot(FingersKinInfo.Finger(3).Es_MaxVel,150,'ob')

figure;
plot(Posq_F3(:,3))
hold on; plot(FingersKinInfo.Finger(3).Fs_MaxPos,-80,'og')
hold on; plot(FingersKinInfo.Finger(3).Es_MaxPos,-50,'or')
hold on; vline(FingersKinInfo.Finger(3).Fs_MaxVel,'g--')
hold on; vline(FingersKinInfo.Finger(3).Es_MaxVel,'r--')

FingersKinData.Finger(3).PurePos=Posq_F3(:,3);
FingersKinData.Finger(3).PureVel=AbsVel_F3;

%% 1-7 F4 movement analysis: Finding the peak points for transition from flextion to extention or extension to flextion

Xq_F4=0:1/(Fs):time_AllTrials{1,4}(end);
Posq_F4=interp1(time_AllTrials{1,4},fingers{1,4}(4).spos,Xq_F4);

half_window=200; %points; about 1 sec: for finding local max and local min in position data

% for finger 1 
[F4_n,F4_m]=size(Xq_F4');

% finding max and min time stamps
j=1;
jj=1;
for i=10*half_window:F4_n-(half_window+1)
    
    logic1=(Posq_F4(i-half_window:i+half_window,3)< Posq_F4(i,3));
    
    if sum(logic1)==2*half_window
        FingersKinInfo.Finger(4).Fs_MaxPos(j,1)=i;
        j=j+1;
        
    end
    
    logic2=(Posq_F4(i-half_window:i+half_window,3)> Posq_F4(i,3));
    
    if sum(logic2)==2*half_window
        FingersKinInfo.Finger(4).Es_MaxPos(jj,1)=i;
        jj=jj+1;
    end
    
    
end

% check on the plot
figure;
plot(Posq_F4(:,3))
hold on; plot(FingersKinInfo.Finger(4).Fs_MaxPos,-80,'og')
hold on; plot(FingersKinInfo.Finger(4).Es_MaxPos,-50,'or')

%% 1-8 F4 movement analysis: Finding the max velocity during transition from flextion to extention or extension to flextion

Velq_F4=interp1(time_AllTrials{1,4},fingers{1,4}(4).svel,Xq_F4);
AbsVel_F4=sqrt(Velq_F4(:,1).^2+Velq_F4(:,2).^2+Velq_F4(:,3).^2);

Index_F4=sort([1; FingersKinInfo.Finger(4).Fs_MaxPos;...
    FingersKinInfo.Finger(4).Es_MaxPos; F4_n]);

j=1;
for i=1:2:(length(Index_F4)-1)
    Section=AbsVel_F4(Index_F4(i):(Index_F4(i+1)));
    [m_S n_S]=max(Section);
    FingersKinInfo.Finger(4).Fs_MaxVel(1,j)=(n_S+Index_F4(i)-1);
    j=j+1;
    
end

j=1;
for i=2:2:(length(Index_F4)-1)
    Section=AbsVel_F4(Index_F4(i):(Index_F4(i+1)));
    [m_S n_S]=max(Section);
    FingersKinInfo.Finger(4).Es_MaxVel(1,j)=(n_S+Index_F4(i)-1);
    j=j+1;
end


% check on the plot
figure;
plot(AbsVel_F4)
hold on; plot(FingersKinInfo.Finger(4).Fs_MaxVel,100,'ok')
hold on; plot(FingersKinInfo.Finger(4).Es_MaxVel,150,'ob')

figure;
plot(Posq_F4(:,3))
hold on; plot(FingersKinInfo.Finger(4).Fs_MaxPos,-80,'og')
hold on; plot(FingersKinInfo.Finger(4).Es_MaxPos,-50,'or')
hold on; vline(FingersKinInfo.Finger(4).Fs_MaxVel,'g--')
hold on; vline(FingersKinInfo.Finger(4).Es_MaxVel,'r--')

FingersKinData.Finger(4).PurePos=Posq_F4(:,3);
FingersKinData.Finger(4).PureVel=AbsVel_F4;

%% 1-9 F5 movement analysis: Finding the peak points for transition from flextion to extention or extension to flextion

Xq_F5=0:1/(Fs):time_AllTrials{1,5}(end);
Posq_F5=interp1(time_AllTrials{1,5},fingers{1,5}(5).spos,Xq_F5);

half_window=350; %points; about 1 sec: for finding local max and local min in position data

% for finger 1 
[F5_n,F5_m]=size(Xq_F5');

% finding max and min time stamps
j=1;
jj=1;
for i=5*half_window:F5_n-(half_window+1)
    
    logic1=(Posq_F5(i-half_window:i+half_window,3)< Posq_F5(i,3));
    
    if sum(logic1)==2*half_window
        FingersKinInfo.Finger(5).Fs_MaxPos(j,1)=i;
        j=j+1;
    end
    
    logic2=(Posq_F5(i-half_window:i+half_window,3)> Posq_F5(i,3));
    
    if sum(logic2)==2*half_window
        FingersKinInfo.Finger(5).Es_MaxPos(jj,1)=i;
        jj=jj+1;
        
    end
     
end

FingersKinInfo.Finger(5).Es_MaxPos=FingersKinInfo.Finger(5).Es_MaxPos([1:16,18:end]);

% check on the plot
figure;
plot(Posq_F5(:,3))
hold on; plot(FingersKinInfo.Finger(5).Fs_MaxPos,-50,'og')
hold on; plot(FingersKinInfo.Finger(5).Es_MaxPos,-20,'or')

%% 1-10 F5 movement analysis: Finding the max velocity during transition from flextion to extention or extension to flextion

Velq_F5=interp1(time_AllTrials{1,5},fingers{1,5}(5).svel,Xq_F5);
AbsVel_F5=sqrt(Velq_F5(:,1).^2+Velq_F5(:,2).^2+Velq_F5(:,3).^2);

Index_F5=sort([1; FingersKinInfo.Finger(5).Fs_MaxPos;...
    FingersKinInfo.Finger(5).Es_MaxPos; F5_n-(half_window+1)]);

j=1;
for i=1:2:(length(Index_F5)-1)
    Section=AbsVel_F5(Index_F5(i):(Index_F5(i+1)));
    [m_S n_S]=max(Section);
    FingersKinInfo.Finger(5).Fs_MaxVel(1,j)=(n_S+Index_F5(i)-1);
    j=j+1;
end

j=1;
for i=2:2:(length(Index_F5)-1)
    Section=AbsVel_F5(Index_F5(i):(Index_F5(i+1)));
    [m_S n_S]=max(Section);
    FingersKinInfo.Finger(5).Es_MaxVel(1,j)=(n_S+Index_F5(i)-1);
    j=j+1;
end

% check on the plot
figure;
plot(AbsVel_F5)
hold on; plot(FingersKinInfo.Finger(5).Fs_MaxVel,100,'ok')
hold on; plot(FingersKinInfo.Finger(5).Es_MaxVel,150,'ob')

figure;
plot(Posq_F5(:,3))
hold on; plot(FingersKinInfo.Finger(5).Fs_MaxPos,-50,'og')
hold on; plot(FingersKinInfo.Finger(5).Es_MaxPos,-20,'or')
hold on; vline(FingersKinInfo.Finger(5).Fs_MaxVel,'g--')
hold on; vline(FingersKinInfo.Finger(5).Es_MaxVel,'r--')

FingersKinData.Finger(5).PurePos=Posq_F5(:,3);
FingersKinData.Finger(5).PureVel=AbsVel_F5;


%% Referencing options for all Fingers

Modification=[length(Xq_F1),length(Xq_F2),length(Xq_F3),length(Xq_F4),length(Xq_F5)];

for Fi=1:5
    % modification to sync
    ECoGdata_Modified=lfp(1:Nch_Record,trial_start(Fi):(trial_start(Fi)+Modification(Fi)-1));
    % z score
    ECoGdata_1=zscore(ECoGdata_Modified')';
    
    % option median or mean:
    switch Reference
        
        case 1
            F_Median =median(ECoGdata_1(Brain_Area,:),1);
            ECoGdata_2= ECoGdata_1-repmat(F_Median,Nch_Record,1);
            
        case 2
            F_Mean =mean(ECoGdata_1(Brain_Area,:),1);
            ECoGdata_2= ECoGdata_1-repmat(F_Mean,Nch_Record,1);
    end
    ECoG_data(Fi).Finger=ECoGdata_2';
     
end

%% save related data
% Fs, FingersKinData, FingersKinInfo, Selected_Chs,    
save('E:\ECoGLeapMotion\ResultsGroupAnalysis\github_Branch_V3/Subject2.mat',...
    'ECoG_data','FingersKinData','FingersKinInfo','Fs','Selected_Chs',...
    'Selected_HandChs', 'Selected_PMChs','Selected_PMSMChs', 'Selected_SMChs');

