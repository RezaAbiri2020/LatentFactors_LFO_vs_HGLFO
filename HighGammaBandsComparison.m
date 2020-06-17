% Analyzing the MLR for all bands and subbands of HG-LFO
% generating the required raw/hilbert signals for LFO and HG-LFO
% save these signals
% in this version; also saving the index for pos max and vel max
% for later sepearting flexion and extension, ERP analysis
% save all these signals/indexes for later phase-coupling
% analysis/SequentialTemporal analysis


close all;
clear all;
clc;


%% General parameters that should be specified at the begining

% choosing brain area to do analysis on that part and do referencing
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

% which band to do analysis: 6 bands
FilterBand=6;

% which feature to do analysis
%Feature=1; % pure filtered bands
Feature=2;   % abs hilbert or envelope


%% loading and breaking raw ECoG data into trials
load('E:\ECoGLeapMotion\DataPatientTwo\ECoGData\ECoG_data.mat');

% Final list of bad channels
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

%% Determining brain section for analysis

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

%% F1 movement analysis: Finding the peak points for transition from flextion to extention or extension to flextion

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

%% F1 movement analysis: Finding the max velocity during transition from flextion to extention or extension to flextion

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


%% F2 movement analysis: Finding the peak points for transition from flextion to extention or extension to flextion

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

%% F2 movement analysis: Finding the max velocity during transition from flextion to extention or extension to flextion

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

%% F3 movement analysis: Finding the peak points for transition from flextion to extention or extension to flextion

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

%% F3 movement analysis: Finding the max velocity during transition from flextion to extention or extension to flextion

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

%% F4 movement analysis: Finding the peak points for transition from flextion to extention or extension to flextion

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

%% F4 movement analysis: Finding the max velocity during transition from flextion to extention or extension to flextion

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

%% F5 movement analysis: Finding the peak points for transition from flextion to extention or extension to flextion

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

%% F5 movement analysis: Finding the max velocity during transition from flextion to extention or extension to flextion

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

%% Referencing & Filtering options for generating raw bands for all Fingers

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
    
    % choose the band freq for analysis
    switch FilterBand
        
        case 1 % delta
            [b,a]=butter(3,[.5, 4]/(Fs/2));
            F_Filtered=filtfilt(b,a,ECoGdata_2');
            Delta_Fingers(Fi).Filtered=F_Filtered;
            F_Hilbert=abs(hilbert(F_Filtered));
            Delta_Fingers(Fi).Hilbert=F_Hilbert;
            
        case 2 % theta
            [b,a]=butter(3,[4, 8]/(Fs/2));
            F_Filtered=filtfilt(b,a,ECoGdata_2');
            F_Hilbert=abs(hilbert(F_Filtered));
            Theta_Fingers(Fi).Hilbert=F_Hilbert;
            
        case 3 % alpha
            [b,a]=butter(3,[8, 13]/(Fs/2));
            F_Filtered=filtfilt(b,a,ECoGdata_2');
            F_Hilbert=abs(hilbert(F_Filtered));
            Alpha_Fingers(Fi).Hilbert=F_Hilbert;
            
        case 4 % beta
            F_FilteredAll=[];
            F_HilbertAll=[];
            
            Betabands={
                [13,19]
                [19,30]};
            
            for band=1:2
                [b,a]=butter(3,Betabands{band}/(Fs/2));
                F_FilteredAll(:,:,band)=filtfilt(b,a,ECoGdata_2');
                
                % Envelope & Hilbert
                [F_Envelope,Lower]=envelope(F_FilteredAll(:,:,band));
                % Hilbert & Filter
                [b,a]=butter(3,[.5, 4]/(Fs/2));
                F_HilbertAll(:,:,band)=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
                
            end
            
            F_Filtered=mean(F_FilteredAll,3);
            F_Hilbert=mean(F_HilbertAll,3);
            Beta_Fingers(Fi).Hilbert=F_Hilbert;
            
        case 5 % low gamma
            LGbands={
                [30,36]
                [36,42]
                [42,50]};
            
            % for whole band
            F_FilteredAll=[];
            F_HilbertAll=[];
            [b,a]=butter(3,[30 50]/(Fs/2));
            F_FilteredAll=filtfilt(b,a,ECoGdata_2');
           
            % Envelope & Hilbert
            [F_Envelope,Lower]=envelope(F_FilteredAll);
            LG_FingersWholeBand(Fi).PureEnv=F_Envelope;
            % Hilbert & Filter in delta
            [b,a]=butter(3,[.5, 4]/(Fs/2));
            F_HilbertAll=filtfilt(b,a,F_Envelope);
            LG_FingersWholeBand(Fi).EnvDeltaWOutMean=F_HilbertAll;
            
            % Hilbert & Filter in theta
            F_FilteredAll=[];
            F_HilbertAll=[];
            [b,a]=butter(3,[4, 8]/(Fs/2));
            F_HilbertAll=filtfilt(b,a,F_Envelope);
            LG_FingersWholeBand(Fi).EnvThetaWOutMean=F_HilbertAll;
            
            % for different sub bands
            AllSubBands=[];
            for band=1:3
                [b,a]=butter(3,LGbands{band}/(Fs/2));
                F_FilteredAll=filtfilt(b,a,ECoGdata_2');
                               
                % Envelope & Hilbert
                [F_Envelope,Lower]=envelope(F_FilteredAll);
                LG_Fingers(Fi).SubBands(band).PureEnv=F_Envelope;
                AllSubBands(:,:,band)=F_Envelope;
                
                % Hilbert & Filter for delta
                [b,a]=butter(3,[.5, 4]/(Fs/2));
                F_HilbertDelta1=filtfilt(b,a,F_Envelope);
                LG_Fingers(Fi).SubBands(band).EnvDeltaWOutMean=F_HilbertDelta1;
               
                % Hilbert & Filter for theta
                [b,a]=butter(3,[4, 8]/(Fs/2));
                F_HilbertTheta1=filtfilt(b,a,F_Envelope);
                LG_Fingers(Fi).SubBands(band).EnvThetaWOutMean=F_HilbertTheta1;
                
            end
            
            % for the averege of bands
            LG_FingersAvgSubBands(Fi).PureEnv=mean(AllSubBands,3);
            % Hilbert & Filter for delta
            [b,a]=butter(3,[.5, 4]/(Fs/2));
            F_HilbertDelta1=filtfilt(b,a,mean(AllSubBands,3));
            LG_FingersAvgSubBands(Fi).EnvDeltaWOutMean=F_HilbertDelta1;
            
            % Hilbert & Filter for theta
            [b,a]=butter(3,[4, 8]/(Fs/2));
            F_HilbertTheta1=filtfilt(b,a,mean(AllSubBands,3));
            LG_FingersAvgSubBands(Fi).EnvThetaWOutMean=F_HilbertTheta1;
            
        case 6
            
            HGbands={
                [70,77]
                [77,85]
                [85,93]
                [93,102]
                [102,113]
                [113,124]
                [124,136]
                [136,150]};
            
            % for whole band
            F_FilteredAll=[];
            F_HilbertAll=[];
            [b,a]=butter(3,[70 150]/(Fs/2));
            F_FilteredAll=filtfilt(b,a,ECoGdata_2');
            HG_FingersWholeBand(Fi).Filtered=F_FilteredAll;
            % Envelope & Hilbert
            [F_Envelope,Lower]=envelope(F_FilteredAll);
            HG_FingersWholeBand(Fi).PureEnv=F_Envelope;
            % Hilbert & Filter in delta
            [b,a]=butter(3,[.5, 4]/(Fs/2));
            F_HilbertAll=filtfilt(b,a,F_Envelope);
            HG_FingersWholeBand(Fi).EnvDeltaWOutMean=F_HilbertAll;
            F_HilbertAll=[];
            F_HilbertAll=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
            HG_FingersWholeBand(Fi).EnvDeltaWithMean=F_HilbertAll;
            
            % Hilbert & Filter in theta
            F_FilteredAll=[];
            F_HilbertAll=[];
            [b,a]=butter(3,[4, 8]/(Fs/2));
            F_HilbertAll=filtfilt(b,a,F_Envelope);
            HG_FingersWholeBand(Fi).EnvThetaWOutMean=F_HilbertAll;
            F_HilbertAll=[];
            F_HilbertAll=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
            HG_FingersWholeBand(Fi).EnvThetaWithMean=F_HilbertAll;
            
            % Hilbert & Filter in alpha
            F_FilteredAll=[];
            F_HilbertAll=[];
            [b,a]=butter(3,[8, 13]/(Fs/2));
            F_HilbertAll=filtfilt(b,a,F_Envelope);
            HG_FingersWholeBand(Fi).EnvAlphaWOutMean=F_HilbertAll;
            F_HilbertAll=[];
            F_HilbertAll=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
            HG_FingersWholeBand(Fi).EnvAlphaWithMean=F_HilbertAll;
            
            % Hilbert & Filter in beta
            F_FilteredAll=[];
            F_HilbertAll=[];
            [b,a]=butter(3,[13, 30]/(Fs/2));
            F_HilbertAll=filtfilt(b,a,F_Envelope);
            HG_FingersWholeBand(Fi).EnvBetaWOutMean=F_HilbertAll;
            F_HilbertAll=[];
            F_HilbertAll=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
            HG_FingersWholeBand(Fi).EnvBetaWithMean=F_HilbertAll;
            
            
            % for different sub bands
            if Fi==1
                for band=1:8
                    [b,a]=butter(3,HGbands{band}/(Fs/2));
                    F_FilteredAll=filtfilt(b,a,ECoGdata_2');
                    HG_Finger1SubBands(band).Filtered=F_FilteredAll;
                    
                    % Envelope & Hilbert
                    [F_Envelope,Lower]=envelope(F_FilteredAll);
                    HG_Finger1SubBands(band).PureEnv=F_Envelope;
                    % Hilbert & Filter for delta
                    [b,a]=butter(3,[.5, 4]/(Fs/2));
                    F_HilbertDelta1=filtfilt(b,a,F_Envelope);
                    F_HilbertDelta2=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
                    HG_Finger1SubBands(band).EnvDeltaWOutMean=F_HilbertDelta1;
                    HG_Finger1SubBands(band).EnvDeltaWithMean=F_HilbertDelta2;
                    
                    % Hilbert & Filter for theta
                    [b,a]=butter(3,[4, 8]/(Fs/2));
                    F_HilbertTheta1=filtfilt(b,a,F_Envelope);
                    F_HilbertTheta2=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
                    HG_Finger1SubBands(band).EnvThetaWOutMean=F_HilbertTheta1;
                    HG_Finger1SubBands(band).EnvThetaWithMean=F_HilbertTheta2;
                    
                    % Hilbert & Filter for alpha
                    [b,a]=butter(3,[8, 13]/(Fs/2));
                    F_HilbertAlpha1=filtfilt(b,a,F_Envelope);
                    F_HilbertAlpha2=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
                    HG_Finger1SubBands(band).EnvAlphaWOutMean=F_HilbertAlpha1;
                    HG_Finger1SubBands(band).EnvAlphaWithMean=F_HilbertAlpha2;
                    
                    % Hilbert & Filter for beta
                    [b,a]=butter(3,[13, 30]/(Fs/2));
                    F_HilbertBeta1=filtfilt(b,a,F_Envelope);
                    F_HilbertBeta2=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
                    HG_Finger1SubBands(band).EnvBetaWOutMean=F_HilbertBeta1;
                    HG_Finger1SubBands(band).EnvBetaWithMean=F_HilbertBeta2;
                    
                end
            elseif Fi==2
                for band=1:8
                    [b,a]=butter(3,HGbands{band}/(Fs/2));
                    F_FilteredAll=filtfilt(b,a,ECoGdata_2');
                    HG_Finger2SubBands(band).Filtered=F_FilteredAll;
                    
                    % Envelope & Hilbert
                    [F_Envelope,Lower]=envelope(F_FilteredAll);
                    HG_Finger2SubBands(band).PureEnv=F_Envelope;
                    % Hilbert & Filter for delta
                    [b,a]=butter(3,[.5, 4]/(Fs/2));
                    F_HilbertDelta1=filtfilt(b,a,F_Envelope);
                    F_HilbertDelta2=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
                    HG_Finger2SubBands(band).EnvDeltaWOutMean=F_HilbertDelta1;
                    HG_Finger2SubBands(band).EnvDeltaWithMean=F_HilbertDelta2;
                    
                    % Hilbert & Filter for theta
                    [b,a]=butter(3,[4, 8]/(Fs/2));
                    F_HilbertTheta1=filtfilt(b,a,F_Envelope);
                    F_HilbertTheta2=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
                    HG_Finger2SubBands(band).EnvThetaWOutMean=F_HilbertTheta1;
                    HG_Finger2SubBands(band).EnvThetaWithMean=F_HilbertTheta2;
                    
                    % Hilbert & Filter for alpha
                    [b,a]=butter(3,[8, 13]/(Fs/2));
                    F_HilbertAlpha1=filtfilt(b,a,F_Envelope);
                    F_HilbertAlpha2=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
                    HG_Finger2SubBands(band).EnvAlphaWOutMean=F_HilbertAlpha1;
                    HG_Finger2SubBands(band).EnvAlphaWithMean=F_HilbertAlpha2;
                    
                    % Hilbert & Filter for beta
                    [b,a]=butter(3,[13, 30]/(Fs/2));
                    F_HilbertBeta1=filtfilt(b,a,F_Envelope);
                    F_HilbertBeta2=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
                    HG_Finger2SubBands(band).EnvBetaWOutMean=F_HilbertBeta1;
                    HG_Finger2SubBands(band).EnvBetaWithMean=F_HilbertBeta2;
                    
                end
            elseif Fi==3
                for band=1:8
                    [b,a]=butter(3,HGbands{band}/(Fs/2));
                    F_FilteredAll=filtfilt(b,a,ECoGdata_2');
                    HG_Finger3SubBands(band).Filtered=F_FilteredAll;
                    
                    % Envelope & Hilbert
                    [F_Envelope,Lower]=envelope(F_FilteredAll);
                    HG_Finger3SubBands(band).PureEnv=F_Envelope;
                    % Hilbert & Filter for delta
                    [b,a]=butter(3,[.5, 4]/(Fs/2));
                    F_HilbertDelta1=filtfilt(b,a,F_Envelope);
                    F_HilbertDelta2=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
                    HG_Finger3SubBands(band).EnvDeltaWOutMean=F_HilbertDelta1;
                    HG_Finger3SubBands(band).EnvDeltaWithMean=F_HilbertDelta2;
                    
                    % Hilbert & Filter for theta
                    [b,a]=butter(3,[4, 8]/(Fs/2));
                    F_HilbertTheta1=filtfilt(b,a,F_Envelope);
                    F_HilbertTheta2=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
                    HG_Finger3SubBands(band).EnvThetaWOutMean=F_HilbertTheta1;
                    HG_Finger3SubBands(band).EnvThetaWithMean=F_HilbertTheta2;
                    
                    % Hilbert & Filter for alpha
                    [b,a]=butter(3,[8, 13]/(Fs/2));
                    F_HilbertAlpha1=filtfilt(b,a,F_Envelope);
                    F_HilbertAlpha2=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
                    HG_Finger3SubBands(band).EnvAlphaWOutMean=F_HilbertAlpha1;
                    HG_Finger3SubBands(band).EnvAlphaWithMean=F_HilbertAlpha2;
                    
                    % Hilbert & Filter for beta
                    [b,a]=butter(3,[13, 30]/(Fs/2));
                    F_HilbertBeta1=filtfilt(b,a,F_Envelope);
                    F_HilbertBeta2=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
                    HG_Finger3SubBands(band).EnvBetaWOutMean=F_HilbertBeta1;
                    HG_Finger3SubBands(band).EnvBetaWithMean=F_HilbertBeta2;
                    
                end
            elseif Fi==4
                for band=1:8
                    [b,a]=butter(3,HGbands{band}/(Fs/2));
                    F_FilteredAll=filtfilt(b,a,ECoGdata_2');
                    HG_Finger4SubBands(band).Filtered=F_FilteredAll;
                    
                    % Envelope & Hilbert
                    [F_Envelope,Lower]=envelope(F_FilteredAll);
                    HG_Finger4SubBands(band).PureEnv=F_Envelope;
                    % Hilbert & Filter for delta
                    [b,a]=butter(3,[.5, 4]/(Fs/2));
                    F_HilbertDelta1=filtfilt(b,a,F_Envelope);
                    F_HilbertDelta2=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
                    HG_Finger4SubBands(band).EnvDeltaWOutMean=F_HilbertDelta1;
                    HG_Finger4SubBands(band).EnvDeltaWithMean=F_HilbertDelta2;
                    
                    % Hilbert & Filter for theta
                    [b,a]=butter(3,[4, 8]/(Fs/2));
                    F_HilbertTheta1=filtfilt(b,a,F_Envelope);
                    F_HilbertTheta2=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
                    HG_Finger4SubBands(band).EnvThetaWOutMean=F_HilbertTheta1;
                    HG_Finger4SubBands(band).EnvThetaWithMean=F_HilbertTheta2;
                    
                    % Hilbert & Filter for alpha
                    [b,a]=butter(3,[8, 13]/(Fs/2));
                    F_HilbertAlpha1=filtfilt(b,a,F_Envelope);
                    F_HilbertAlpha2=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
                    HG_Finger4SubBands(band).EnvAlphaWOutMean=F_HilbertAlpha1;
                    HG_Finger4SubBands(band).EnvAlphaWithMean=F_HilbertAlpha2;
                    
                    % Hilbert & Filter for beta
                    [b,a]=butter(3,[13, 30]/(Fs/2));
                    F_HilbertBeta1=filtfilt(b,a,F_Envelope);
                    F_HilbertBeta2=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
                    HG_Finger4SubBands(band).EnvBetaWOutMean=F_HilbertBeta1;
                    HG_Finger4SubBands(band).EnvBetaWithMean=F_HilbertBeta2;
                    
                end
            elseif Fi==5
                for band=1:8
                    [b,a]=butter(3,HGbands{band}/(Fs/2));
                    F_FilteredAll=filtfilt(b,a,ECoGdata_2');
                    HG_Finger5SubBands(band).Filtered=F_FilteredAll;
                    
                    % Envelope & Hilbert
                    [F_Envelope,Lower]=envelope(F_FilteredAll);
                    HG_Finger5SubBands(band).PureEnv=F_Envelope;
                    % Hilbert & Filter for delta
                    [b,a]=butter(3,[.5, 4]/(Fs/2));
                    F_HilbertDelta1=filtfilt(b,a,F_Envelope);
                    F_HilbertDelta2=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
                    HG_Finger5SubBands(band).EnvDeltaWOutMean=F_HilbertDelta1;
                    HG_Finger5SubBands(band).EnvDeltaWithMean=F_HilbertDelta2;
                    
                    % Hilbert & Filter for theta
                    [b,a]=butter(3,[4, 8]/(Fs/2));
                    F_HilbertTheta1=filtfilt(b,a,F_Envelope);
                    F_HilbertTheta2=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
                    HG_Finger5SubBands(band).EnvThetaWOutMean=F_HilbertTheta1;
                    HG_Finger5SubBands(band).EnvThetaWithMean=F_HilbertTheta2;
                    
                    % Hilbert & Filter for alpha
                    [b,a]=butter(3,[8, 13]/(Fs/2));
                    F_HilbertAlpha1=filtfilt(b,a,F_Envelope);
                    F_HilbertAlpha2=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
                    HG_Finger5SubBands(band).EnvAlphaWOutMean=F_HilbertAlpha1;
                    HG_Finger5SubBands(band).EnvAlphaWithMean=F_HilbertAlpha2;
                    
                    % Hilbert & Filter for beta
                    [b,a]=butter(3,[13, 30]/(Fs/2));
                    F_HilbertBeta1=filtfilt(b,a,F_Envelope);
                    F_HilbertBeta2=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
                    HG_Finger5SubBands(band).EnvBetaWOutMean=F_HilbertBeta1;
                    HG_Finger5SubBands(band).EnvBetaWithMean=F_HilbertBeta2;
                    
                end
            end      
    end
      
     
end

%save('E:\ECoGLeapMotion\DataPatientTwo\github_Branch_V1/LFO_signals.mat','LFO_signals')
save('E:\ECoGLeapMotion\DataPatientTwo\github_Branch_V1/FingersKinIndexes.mat','FingersKinIndexes')
%save('E:\ECoGLeapMotion\DataPatientTwo\github_Branch_V1/HG_Direct_LFO_Signals.mat','HG_Direct_LFO_Signals')

%% if HG Continued: similar analysis for the avereged subbands of highgamma 
Fi=1;
AllSubBands=[];
for band=1:8
    AllSubBands(:,:,band)=HG_Finger1SubBands(band).PureEnv;
end
HG_FingersAvgSubBands(Fi).PureEnv=mean(AllSubBands,3);

Fi=2;
AllSubBands=[];
for band=1:8
    AllSubBands(:,:,band)=HG_Finger2SubBands(band).PureEnv;
end
HG_FingersAvgSubBands(Fi).PureEnv=mean(AllSubBands,3);

Fi=3;
AllSubBands=[];
for band=1:8
    AllSubBands(:,:,band)=HG_Finger3SubBands(band).PureEnv;
end
HG_FingersAvgSubBands(Fi).PureEnv=mean(AllSubBands,3);

Fi=4;
AllSubBands=[];
for band=1:8
    AllSubBands(:,:,band)=HG_Finger4SubBands(band).PureEnv;
end
HG_FingersAvgSubBands(Fi).PureEnv=mean(AllSubBands,3);

Fi=5;
AllSubBands=[];
for band=1:8
    AllSubBands(:,:,band)=HG_Finger5SubBands(band).PureEnv;
end
HG_FingersAvgSubBands(Fi).PureEnv=mean(AllSubBands,3);

for Fi=1:5
      
    % Hilbert & Filter in delta
    PureEnv=[];
    PureEnv=HG_FingersAvgSubBands(Fi).PureEnv;
    [b,a]=butter(3,[.5, 4]/(Fs/2));
    HG_FingersAvgSubBands(Fi).EnvDeltaWOutMean=filtfilt(b,a,PureEnv);
    HG_FingersAvgSubBands(Fi).EnvDeltaWithMean=...
        filtfilt(b,a,PureEnv)+repmat(mean(PureEnv),size(PureEnv,1),1);
    
    % Hilbert & Filter in Theta
    [b,a]=butter(3,[4, 8]/(Fs/2));
    HG_FingersAvgSubBands(Fi).EnvThetaWOutMean=filtfilt(b,a,PureEnv);
    HG_FingersAvgSubBands(Fi).EnvThetaWithMean=...
        filtfilt(b,a,PureEnv)+repmat(mean(PureEnv),size(PureEnv,1),1);
    
    % Hilbert & Filter in Alpha
    [b,a]=butter(3,[8, 13]/(Fs/2));
    HG_FingersAvgSubBands(Fi).EnvAlphaWOutMean=filtfilt(b,a,PureEnv);
    HG_FingersAvgSubBands(Fi).EnvAlphaWithMean=...
        filtfilt(b,a,PureEnv)+repmat(mean(PureEnv),size(PureEnv,1),1);
    
    % Hilbert & Filter in Beta
    [b,a]=butter(3,[13, 30]/(Fs/2));
    HG_FingersAvgSubBands(Fi).EnvBetaWOutMean=filtfilt(b,a,PureEnv);
    HG_FingersAvgSubBands(Fi).EnvBetaWithMean=...
        filtfilt(b,a,PureEnv)+repmat(mean(PureEnv),size(PureEnv,1),1);
    
end

%save('E:\ECoGLeapMotion\DataPatientTwo\github_Branch_V1/HG_Avg_LFO_Signals.mat','HG_Avg_LFO_Signals')

%% Multiple Linear Regression for delta band
    
Kin_Pos={Posq_F1(:,1) Posq_F2(:,3) Posq_F3(:,3) Posq_F4(:,3) Posq_F5(:,3)};
Kin_Vel={AbsVel_F1 AbsVel_F2 AbsVel_F3 AbsVel_F4 AbsVel_F5};

% for delta band
% Input_dynamics=Delta_Fingers(Fi).Hilbert=F_Hilbert
% show the plots?
Figures=1;

% prediction the velocities and positions and cacculating the R2
R2FingersVel=[];
R2FingersPos=[];
for Fi=1:5
    Input_dynamics=Delta_Fingers(Fi).Hilbert;
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2FingersVel=[R2FingersVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2FingersPos=[R2FingersPos,R2Pos];
end
Delta_FingersR2Vel=R2FingersVel;
Delta_FingersR2Pos=R2FingersPos;

%% save R2s for delta

save('Delta_AllR2s.mat','Delta_FingersR2Vel','Delta_FingersR2Pos');

%% Multiple Linear Regression for theta band
    
Kin_Pos={Posq_F1(:,1) Posq_F2(:,3) Posq_F3(:,3) Posq_F4(:,3) Posq_F5(:,3)};
Kin_Vel={AbsVel_F1 AbsVel_F2 AbsVel_F3 AbsVel_F4 AbsVel_F5};

% for delta band
% Input_dynamics=Delta_Fingers(Fi).Hilbert=F_Hilbert
% show the plots?
Figures=1;

% prediction the velocities and positions and cacculating the R2
R2FingersVel=[];
R2FingersPos=[];
for Fi=1:5
    Input_dynamics=Theta_Fingers(Fi).Hilbert;
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2FingersVel=[R2FingersVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2FingersPos=[R2FingersPos,R2Pos];
end
Theta_FingersR2Vel=R2FingersVel;
Theta_FingersR2Pos=R2FingersPos;

%% save R2s for theta

save('theta_AllR2s.mat','Theta_FingersR2Vel','Theta_FingersR2Pos');

%% Multiple Linear Regression for alpha band
    
Kin_Pos={Posq_F1(:,1) Posq_F2(:,3) Posq_F3(:,3) Posq_F4(:,3) Posq_F5(:,3)};
Kin_Vel={AbsVel_F1 AbsVel_F2 AbsVel_F3 AbsVel_F4 AbsVel_F5};

% for delta band
% Input_dynamics=Delta_Fingers(Fi).Hilbert=F_Hilbert
% show the plots?
Figures=1;

% prediction the velocities and positions and cacculating the R2
R2FingersVel=[];
R2FingersPos=[];
for Fi=1:5
    Input_dynamics=Alpha_Fingers(Fi).Hilbert;
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2FingersVel=[R2FingersVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2FingersPos=[R2FingersPos,R2Pos];
end
Alpha_FingersR2Vel=R2FingersVel;
Alpha_FingersR2Pos=R2FingersPos;

%% save R2s for alpha

save('Alpha_AllR2s.mat','Alpha_FingersR2Vel','Alpha_FingersR2Pos');

%% Multiple Linear Regression for beta band
    
Kin_Pos={Posq_F1(:,1) Posq_F2(:,3) Posq_F3(:,3) Posq_F4(:,3) Posq_F5(:,3)};
Kin_Vel={AbsVel_F1 AbsVel_F2 AbsVel_F3 AbsVel_F4 AbsVel_F5};

% for delta band
% Input_dynamics=Delta_Fingers(Fi).Hilbert=F_Hilbert
% show the plots?
Figures=1;

% prediction the velocities and positions and cacculating the R2
R2FingersVel=[];
R2FingersPos=[];
for Fi=1:5
    Input_dynamics=Beta_Fingers(Fi).Hilbert;
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2FingersVel=[R2FingersVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2FingersPos=[R2FingersPos,R2Pos];
end
Beta_FingersR2Vel=R2FingersVel;
Beta_FingersR2Pos=R2FingersPos;

%% save R2s for beta

save('Beta_AllR2s.mat','Beta_FingersR2Vel','Beta_FingersR2Pos');

%% Multiple Linear Regression for low gamma band
    
Kin_Pos={Posq_F1(:,1) Posq_F2(:,3) Posq_F3(:,3) Posq_F4(:,3) Posq_F5(:,3)};
Kin_Vel={AbsVel_F1 AbsVel_F2 AbsVel_F3 AbsVel_F4 AbsVel_F5};

% show the plots?
Figures=1;

% prediction the velocities and positions and cacculating the R2

% for whole band
R2FingersVel=[];
R2FingersPos=[];
for Fi=1:5
    Input_dynamics=LG_FingersWholeBand(Fi).PureEnv;
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2FingersVel=[R2FingersVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2FingersPos=[R2FingersPos,R2Pos];
end
LG_FingersWholeBandR2Vel.PureEnv=R2FingersVel;
LG_FingersWholeBandR2Pos.PureEnv=R2FingersPos;

R2FingersVel=[];
R2FingersPos=[];
for Fi=1:5
    Input_dynamics=LG_FingersWholeBand(Fi).EnvDeltaWOutMean;
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2FingersVel=[R2FingersVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2FingersPos=[R2FingersPos,R2Pos];
end
LG_FingersWholeBandR2Vel.EnvDeltaWOutMean=R2FingersVel;
LG_FingersWholeBandR2Pos.EnvDeltaWOutMean=R2FingersPos;

R2FingersVel=[];
R2FingersPos=[];
for Fi=1:5
    Input_dynamics=LG_FingersWholeBand(Fi).EnvThetaWOutMean;
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2FingersVel=[R2FingersVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2FingersPos=[R2FingersPos,R2Pos];
end
LG_FingersWholeBandR2Vel.EnvThetaWOutMean=R2FingersVel;
LG_FingersWholeBandR2Pos.EnvThetaWOutMean=R2FingersPos;

% for subbands
for Fi=1:5
    R2SubBandsVel=[];
    R2SubBandsPos=[];
    for band=1:3
        Input_dynamics=LG_Fingers(Fi).SubBands(band).PureEnv;
        Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
        Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
        R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
        R2SubBandsVel=[R2SubBandsVel,R2Vel];
        R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
        R2SubBandsPos=[R2SubBandsPos,R2Pos];
    end
    LG_FingersSubBandsR2Vel(Fi).PureEnv=R2SubBandsVel;
    LG_FingersSubBandsR2Pos(Fi).PureEnv=R2SubBandsPos;
end

for Fi=1:5
    R2SubBandsVel=[];
    R2SubBandsPos=[];
    for band=1:3
        Input_dynamics=LG_Fingers(Fi).SubBands(band).EnvDeltaWOutMean;
        Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
        Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
        R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
        R2SubBandsVel=[R2SubBandsVel,R2Vel];
        R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
        R2SubBandsPos=[R2SubBandsPos,R2Pos];
    end
    LG_FingersSubBandsR2Vel(Fi).EnvDeltaWOutMean=R2SubBandsVel;
    LG_FingersSubBandsR2Pos(Fi).EnvDeltaWOutMean=R2SubBandsPos;
end

for Fi=1:5
    R2SubBandsVel=[];
    R2SubBandsPos=[];
    for band=1:3
        Input_dynamics=LG_Fingers(Fi).SubBands(band).EnvThetaWOutMean;
        Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
        Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
        R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
        R2SubBandsVel=[R2SubBandsVel,R2Vel];
        R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
        R2SubBandsPos=[R2SubBandsPos,R2Pos];
    end
    LG_FingersSubBandsR2Vel(Fi).EnvThetaWOutMean=R2SubBandsVel;
    LG_FingersSubBandsR2Pos(Fi).EnvThetaWOutMean=R2SubBandsPos;
end

% for the avg band
R2FingersVel=[];
R2FingersPos=[];
for Fi=1:5
    Input_dynamics=LG_FingersAvgSubBands(Fi).PureEnv;
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2FingersVel=[R2FingersVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2FingersPos=[R2FingersPos,R2Pos];
end
LG_FingersAvgSubBandsR2Vel.PureEnv=R2FingersVel;
LG_FingersAvgSubBandsR2Pos.PureEnv=R2FingersPos;

R2FingersVel=[];
R2FingersPos=[];
for Fi=1:5
    Input_dynamics=LG_FingersAvgSubBands(Fi).EnvDeltaWOutMean;
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2FingersVel=[R2FingersVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2FingersPos=[R2FingersPos,R2Pos];
end
LG_FingersAvgSubBandsR2Vel.EnvDeltaWOutMean=R2FingersVel;
LG_FingersAvgSubBandsR2Pos.EnvDeltaWOutMean=R2FingersPos;

R2FingersVel=[];
R2FingersPos=[];
for Fi=1:5
    Input_dynamics=LG_FingersAvgSubBands(Fi).EnvThetaWOutMean;
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2FingersVel=[R2FingersVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2FingersPos=[R2FingersPos,R2Pos];
end
LG_FingersAvgSubBandsR2Vel.EnvThetaWOutMean=R2FingersVel;
LG_FingersAvgSubBandsR2Pos.EnvThetaWOutMean=R2FingersPos;



%% save R2s for low gamma band 

save('LG_AllR2s.mat','LG_FingersWholeBandR2Vel','LG_FingersWholeBandR2Pos',...
    'LG_FingersSubBandsR2Vel','LG_FingersSubBandsR2Pos',...
    'LG_FingersAvgSubBandsR2Vel','LG_FingersAvgSubBandsR2Pos');

%% Multiple Linear Regression for whole high gamma band
Kin_Pos={Posq_F1(:,1) Posq_F2(:,3) Posq_F3(:,3) Posq_F4(:,3) Posq_F5(:,3)};
Kin_Vel={AbsVel_F1 AbsVel_F2 AbsVel_F3 AbsVel_F4 AbsVel_F5};

% show the plots?
Figures=1;

% Input_dynamics=HG_FingersWholeBand(Fi)....
% prediction the velocities and positions and cacculating the R2
R2FingersVel=[];
R2FingersPos=[];
for Fi=1:5
    Input_dynamics=HG_FingersWholeBand(Fi).Filtered;
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2FingersVel=[R2FingersVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2FingersPos=[R2FingersPos,R2Pos];
end
HG_FingersWholeBandR2Vel.Filtered=R2FingersVel;
HG_FingersWholeBandR2Pos.Filtered=R2FingersPos;

R2FingersVel=[];
R2FingersPos=[];
for Fi=1:5
    Input_dynamics=HG_FingersWholeBand(Fi).PureEnv;
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2FingersVel=[R2FingersVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2FingersPos=[R2FingersPos,R2Pos];
end
HG_FingersWholeBandR2Vel.PureEnv=R2FingersVel;
HG_FingersWholeBandR2Pos.PureEnv=R2FingersPos;

R2FingersVel=[];
R2FingersPos=[];
for Fi=1:5
    Input_dynamics=abs(hilbert(HG_FingersWholeBand(Fi).EnvDeltaWOutMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2FingersVel=[R2FingersVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2FingersPos=[R2FingersPos,R2Pos];
end
HG_FingersWholeBandR2Vel.EnvDeltaWOutMean=R2FingersVel;
HG_FingersWholeBandR2Pos.EnvDeltaWOutMean=R2FingersPos;

R2FingersVel=[];
R2FingersPos=[];
for Fi=1:5
    Input_dynamics=abs(hilbert(HG_FingersWholeBand(Fi).EnvDeltaWithMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2FingersVel=[R2FingersVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2FingersPos=[R2FingersPos,R2Pos];
end
HG_FingersWholeBandR2Vel.EnvDeltaWithMean=R2FingersVel;
HG_FingersWholeBandR2Pos.EnvDeltaWithMean=R2FingersPos;

R2FingersVel=[];
R2FingersPos=[];
for Fi=1:5
    Input_dynamics=abs(hilbert(HG_FingersWholeBand(Fi).EnvThetaWOutMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2FingersVel=[R2FingersVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2FingersPos=[R2FingersPos,R2Pos];
end
HG_FingersWholeBandR2Vel.EnvThetaWOutMean=R2FingersVel;
HG_FingersWholeBandR2Pos.EnvThetaWOutMean=R2FingersPos;

R2FingersVel=[];
R2FingersPos=[];
for Fi=1:5
    Input_dynamics=abs(hilbert(HG_FingersWholeBand(Fi).EnvThetaWithMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2FingersVel=[R2FingersVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2FingersPos=[R2FingersPos,R2Pos];
end
HG_FingersWholeBandR2Vel.EnvThetaWithMean=R2FingersVel;
HG_FingersWholeBandR2Pos.EnvThetaWithMean=R2FingersPos;

R2FingersVel=[];
R2FingersPos=[];
for Fi=1:5
    Input_dynamics=abs(hilbert(HG_FingersWholeBand(Fi).EnvAlphaWOutMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2FingersVel=[R2FingersVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2FingersPos=[R2FingersPos,R2Pos];
end
HG_FingersWholeBandR2Vel.EnvAlphaWOutMean=R2FingersVel;
HG_FingersWholeBandR2Pos.EnvAlphaWOutMean=R2FingersPos;

R2FingersVel=[];
R2FingersPos=[];
for Fi=1:5
    Input_dynamics=abs(hilbert(HG_FingersWholeBand(Fi).EnvAlphaWithMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2FingersVel=[R2FingersVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2FingersPos=[R2FingersPos,R2Pos];
end
HG_FingersWholeBandR2Vel.EnvAlphaWithMean=R2FingersVel;
HG_FingersWholeBandR2Pos.EnvAlphaWithMean=R2FingersPos;

R2FingersVel=[];
R2FingersPos=[];
for Fi=1:5
    Input_dynamics=abs(hilbert(HG_FingersWholeBand(Fi).EnvBetaWOutMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2FingersVel=[R2FingersVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2FingersPos=[R2FingersPos,R2Pos];
end
HG_FingersWholeBandR2Vel.EnvBetaWOutMean=R2FingersVel;
HG_FingersWholeBandR2Pos.EnvBetaWOutMean=R2FingersPos;

R2FingersVel=[];
R2FingersPos=[];
for Fi=1:5
    Input_dynamics=abs(hilbert(HG_FingersWholeBand(Fi).EnvBetaWithMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2FingersVel=[R2FingersVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2FingersPos=[R2FingersPos,R2Pos];
end
HG_FingersWholeBandR2Vel.EnvBetaWithMean=R2FingersVel;
HG_FingersWholeBandR2Pos.EnvBetaWithMean=R2FingersPos;


%% Multiple Linear Regression for individual subbands of high gamma
% show the plots?
Figures=1;
% Input_dynamics=HG_Finger?SubBands....
% prediction the velocities and positions and calculating R2

%%%%%%%%%%%%%%%%%%%%%%% Finger 1
Fi=1;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=HG_Finger1SubBands(band).Filtered;
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).Filtered=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).Filtered=R2SubBandsPos;

Fi=1;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=HG_Finger1SubBands(band).PureEnv;
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).PureEnv=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).PureEnv=R2SubBandsPos;

Fi=1;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=abs(hilbert(HG_Finger1SubBands(band).EnvDeltaWOutMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).EnvDeltaWOutMean=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).EnvDeltaWOutMean=R2SubBandsPos;

Fi=1;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=abs(hilbert(HG_Finger1SubBands(band).EnvDeltaWithMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).EnvDeltaWithMean=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).EnvDeltaWithMean=R2SubBandsPos;

Fi=1;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=abs(hilbert(HG_Finger1SubBands(band).EnvThetaWOutMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).EnvThetaWOutMean=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).EnvThetaWOutMean=R2SubBandsPos;

Fi=1;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=abs(hilbert(HG_Finger1SubBands(band).EnvThetaWithMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).EnvThetaWithMean=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).EnvThetaWithMean=R2SubBandsPos;

Fi=1;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=abs(hilbert(HG_Finger1SubBands(band).EnvAlphaWOutMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).EnvAlphaWOutMean=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).EnvAlphaWOutMean=R2SubBandsPos;

Fi=1;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=abs(hilbert(HG_Finger1SubBands(band).EnvAlphaWithMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).EnvAlphaWithMean=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).EnvAlphaWithMean=R2SubBandsPos;

Fi=1;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=abs(hilbert(HG_Finger1SubBands(band).EnvBetaWOutMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).EnvBetaWOutMean=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).EnvBetaWOutMean=R2SubBandsPos;

Fi=1;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=abs(hilbert(HG_Finger1SubBands(band).EnvBetaWithMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).EnvBetaWithMean=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).EnvBetaWithMean=R2SubBandsPos;

%%%%%%%%%%%%%%%%%%%%%%% Finger 2
Fi=2;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=HG_Finger2SubBands(band).Filtered;
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).Filtered=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).Filtered=R2SubBandsPos;

Fi=2;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=HG_Finger2SubBands(band).PureEnv;
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).PureEnv=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).PureEnv=R2SubBandsPos;

Fi=2;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=abs(hilbert(HG_Finger2SubBands(band).EnvDeltaWOutMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).EnvDeltaWOutMean=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).EnvDeltaWOutMean=R2SubBandsPos;

Fi=2;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=abs(hilbert(HG_Finger2SubBands(band).EnvDeltaWithMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).EnvDeltaWithMean=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).EnvDeltaWithMean=R2SubBandsPos;

Fi=2;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=abs(hilbert(HG_Finger2SubBands(band).EnvThetaWOutMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).EnvThetaWOutMean=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).EnvThetaWOutMean=R2SubBandsPos;

Fi=2;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=abs(hilbert(HG_Finger2SubBands(band).EnvThetaWithMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).EnvThetaWithMean=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).EnvThetaWithMean=R2SubBandsPos;

Fi=2;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=abs(hilbert(HG_Finger2SubBands(band).EnvAlphaWOutMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).EnvAlphaWOutMean=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).EnvAlphaWOutMean=R2SubBandsPos;

Fi=2;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=abs(hilbert(HG_Finger2SubBands(band).EnvAlphaWithMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).EnvAlphaWithMean=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).EnvAlphaWithMean=R2SubBandsPos;

Fi=2;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=abs(hilbert(HG_Finger2SubBands(band).EnvBetaWOutMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).EnvBetaWOutMean=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).EnvBetaWOutMean=R2SubBandsPos;

Fi=2;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=abs(hilbert(HG_Finger2SubBands(band).EnvBetaWithMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).EnvBetaWithMean=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).EnvBetaWithMean=R2SubBandsPos;

%%%%%%%%%%%%%%%%%%%%%%% Finger 3
Fi=3;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=HG_Finger3SubBands(band).Filtered;
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).Filtered=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).Filtered=R2SubBandsPos;

Fi=3;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=HG_Finger3SubBands(band).PureEnv;
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).PureEnv=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).PureEnv=R2SubBandsPos;

Fi=3;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=abs(hilbert(HG_Finger3SubBands(band).EnvDeltaWOutMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).EnvDeltaWOutMean=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).EnvDeltaWOutMean=R2SubBandsPos;

Fi=3;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=abs(hilbert(HG_Finger3SubBands(band).EnvDeltaWithMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).EnvDeltaWithMean=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).EnvDeltaWithMean=R2SubBandsPos;

Fi=3;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=abs(hilbert(HG_Finger3SubBands(band).EnvThetaWOutMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).EnvThetaWOutMean=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).EnvThetaWOutMean=R2SubBandsPos;

Fi=3;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=abs(hilbert(HG_Finger3SubBands(band).EnvThetaWithMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).EnvThetaWithMean=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).EnvThetaWithMean=R2SubBandsPos;

Fi=3;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=abs(hilbert(HG_Finger3SubBands(band).EnvAlphaWOutMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).EnvAlphaWOutMean=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).EnvAlphaWOutMean=R2SubBandsPos;

Fi=3;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=abs(hilbert(HG_Finger3SubBands(band).EnvAlphaWithMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).EnvAlphaWithMean=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).EnvAlphaWithMean=R2SubBandsPos;

Fi=3;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=abs(hilbert(HG_Finger3SubBands(band).EnvBetaWOutMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).EnvBetaWOutMean=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).EnvBetaWOutMean=R2SubBandsPos;

Fi=3;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=abs(hilbert(HG_Finger3SubBands(band).EnvBetaWithMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).EnvBetaWithMean=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).EnvBetaWithMean=R2SubBandsPos;

%%%%%%%%%%%%%%%%%%%%%%% Finger 4
Fi=4;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=HG_Finger4SubBands(band).Filtered;
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).Filtered=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).Filtered=R2SubBandsPos;

Fi=4;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=HG_Finger4SubBands(band).PureEnv;
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).PureEnv=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).PureEnv=R2SubBandsPos;

Fi=4;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=abs(hilbert(HG_Finger4SubBands(band).EnvDeltaWOutMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).EnvDeltaWOutMean=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).EnvDeltaWOutMean=R2SubBandsPos;

Fi=4;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=abs(hilbert(HG_Finger4SubBands(band).EnvDeltaWithMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).EnvDeltaWithMean=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).EnvDeltaWithMean=R2SubBandsPos;

Fi=4;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=abs(hilbert(HG_Finger4SubBands(band).EnvThetaWOutMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).EnvThetaWOutMean=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).EnvThetaWOutMean=R2SubBandsPos;

Fi=4;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=abs(hilbert(HG_Finger4SubBands(band).EnvThetaWithMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).EnvThetaWithMean=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).EnvThetaWithMean=R2SubBandsPos;

Fi=4;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=abs(hilbert(HG_Finger4SubBands(band).EnvAlphaWOutMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).EnvAlphaWOutMean=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).EnvAlphaWOutMean=R2SubBandsPos;

Fi=4;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=abs(hilbert(HG_Finger4SubBands(band).EnvAlphaWithMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).EnvAlphaWithMean=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).EnvAlphaWithMean=R2SubBandsPos;

Fi=4;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=abs(hilbert(HG_Finger4SubBands(band).EnvBetaWOutMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).EnvBetaWOutMean=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).EnvBetaWOutMean=R2SubBandsPos;

Fi=4;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=abs(hilbert(HG_Finger4SubBands(band).EnvBetaWithMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).EnvBetaWithMean=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).EnvBetaWithMean=R2SubBandsPos;

%%%%%%%%%%%%%%%%%%%%%%% Finger 5
Fi=5;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=HG_Finger5SubBands(band).Filtered;
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).Filtered=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).Filtered=R2SubBandsPos;

Fi=5;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=HG_Finger5SubBands(band).PureEnv;
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).PureEnv=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).PureEnv=R2SubBandsPos;

Fi=5;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=abs(hilbert(HG_Finger5SubBands(band).EnvDeltaWOutMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).EnvDeltaWOutMean=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).EnvDeltaWOutMean=R2SubBandsPos;

Fi=5;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=abs(hilbert(HG_Finger5SubBands(band).EnvDeltaWithMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).EnvDeltaWithMean=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).EnvDeltaWithMean=R2SubBandsPos;

Fi=5;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=abs(hilbert(HG_Finger5SubBands(band).EnvThetaWOutMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).EnvThetaWOutMean=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).EnvThetaWOutMean=R2SubBandsPos;

Fi=5;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=abs(hilbert(HG_Finger5SubBands(band).EnvThetaWithMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).EnvThetaWithMean=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).EnvThetaWithMean=R2SubBandsPos;

Fi=5;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=abs(hilbert(HG_Finger5SubBands(band).EnvAlphaWOutMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).EnvAlphaWOutMean=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).EnvAlphaWOutMean=R2SubBandsPos;

Fi=5;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=abs(hilbert(HG_Finger5SubBands(band).EnvAlphaWithMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).EnvAlphaWithMean=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).EnvAlphaWithMean=R2SubBandsPos;

Fi=5;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=abs(hilbert(HG_Finger5SubBands(band).EnvBetaWOutMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).EnvBetaWOutMean=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).EnvBetaWOutMean=R2SubBandsPos;

Fi=5;
R2SubBandsVel=[];
R2SubBandsPos=[];
for band=1:8
    Input_dynamics=abs(hilbert(HG_Finger5SubBands(band).EnvBetaWithMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2SubBandsVel=[R2SubBandsVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2SubBandsPos=[R2SubBandsPos,R2Pos];
end
HG_FingersSubBandsR2Vel(Fi).EnvBetaWithMean=R2SubBandsVel;
HG_FingersSubBandsR2Pos(Fi).EnvBetaWithMean=R2SubBandsPos;

%% Multiple Linear Regression for for the average of subbands of high gamma
% 
% Input_dynamics=HG_FingersAvgSubBands(?)....
% show the plots?
Figures=1;

% Input_dynamics=HG_FingersAvgSubBands(Fi)....
% prediction the velocities and positions and cacculating the R2
R2FingersVel=[];
R2FingersPos=[];
for Fi=1:5
    Input_dynamics=HG_FingersAvgSubBands(Fi).PureEnv;
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2FingersVel=[R2FingersVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2FingersPos=[R2FingersPos,R2Pos];
end
HG_FingersAvgSubBandsR2Vel.PureEnv=R2FingersVel;
HG_FingersAvgSubBandsR2Pos.PureEnv=R2FingersPos;

R2FingersVel=[];
R2FingersPos=[];
for Fi=1:5
    Input_dynamics=abs(hilbert(HG_FingersAvgSubBands(Fi).EnvDeltaWOutMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2FingersVel=[R2FingersVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2FingersPos=[R2FingersPos,R2Pos];
end
HG_FingersAvgSubBandsR2Vel.EnvDeltaWOutMean=R2FingersVel;
HG_FingersAvgSubBandsR2Pos.EnvDeltaWOutMean=R2FingersPos;

R2FingersVel=[];
R2FingersPos=[];
for Fi=1:5
    Input_dynamics=abs(hilbert(HG_FingersAvgSubBands(Fi).EnvDeltaWithMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2FingersVel=[R2FingersVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2FingersPos=[R2FingersPos,R2Pos];
end
HG_FingersAvgSubBandsR2Vel.EnvDeltaWithMean=R2FingersVel;
HG_FingersAvgSubBandsR2Pos.EnvDeltaWithMean=R2FingersPos;

R2FingersVel=[];
R2FingersPos=[];
for Fi=1:5
    Input_dynamics=abs(hilbert(HG_FingersAvgSubBands(Fi).EnvThetaWOutMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2FingersVel=[R2FingersVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2FingersPos=[R2FingersPos,R2Pos];
end
HG_FingersAvgSubBandsR2Vel.EnvThetaWOutMean=R2FingersVel;
HG_FingersAvgSubBandsR2Pos.EnvThetaWOutMean=R2FingersPos;

R2FingersVel=[];
R2FingersPos=[];
for Fi=1:5
    Input_dynamics=abs(hilbert(HG_FingersAvgSubBands(Fi).EnvThetaWithMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2FingersVel=[R2FingersVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2FingersPos=[R2FingersPos,R2Pos];
end
HG_FingersAvgSubBandsR2Vel.EnvThetaWithMean=R2FingersVel;
HG_FingersAvgSubBandsR2Pos.EnvThetaWithMean=R2FingersPos;

R2FingersVel=[];
R2FingersPos=[];
for Fi=1:5
    Input_dynamics=abs(hilbert(HG_FingersAvgSubBands(Fi).EnvAlphaWOutMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2FingersVel=[R2FingersVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2FingersPos=[R2FingersPos,R2Pos];
end
HG_FingersAvgSubBandsR2Vel.EnvAlphaWOutMean=R2FingersVel;
HG_FingersAvgSubBandsR2Pos.EnvAlphaWOutMean=R2FingersPos;

R2FingersVel=[];
R2FingersPos=[];
for Fi=1:5
    Input_dynamics=abs(hilbert(HG_FingersAvgSubBands(Fi).EnvAlphaWithMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2FingersVel=[R2FingersVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2FingersPos=[R2FingersPos,R2Pos];
end
HG_FingersAvgSubBandsR2Vel.EnvAlphaWithMean=R2FingersVel;
HG_FingersAvgSubBandsR2Pos.EnvAlphaWithMean=R2FingersPos;

R2FingersVel=[];
R2FingersPos=[];
for Fi=1:5
    Input_dynamics=abs(hilbert(HG_FingersAvgSubBands(Fi).EnvBetaWOutMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2FingersVel=[R2FingersVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2FingersPos=[R2FingersPos,R2Pos];
end
HG_FingersAvgSubBandsR2Vel.EnvBetaWOutMean=R2FingersVel;
HG_FingersAvgSubBandsR2Pos.EnvBetaWOutMean=R2FingersPos;

R2FingersVel=[];
R2FingersPos=[];
for Fi=1:5
    Input_dynamics=abs(hilbert(HG_FingersAvgSubBands(Fi).EnvBetaWithMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    R2Vel=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    R2FingersVel=[R2FingersVel,R2Vel];
    R2Pos=MultipleRegFunc(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    R2FingersPos=[R2FingersPos,R2Pos];
end
HG_FingersAvgSubBandsR2Vel.EnvBetaWithMean=R2FingersVel;
HG_FingersAvgSubBandsR2Pos.EnvBetaWithMean=R2FingersPos;

%% save R2s for HG

save('HG_AllR2s.mat','HG_FingersWholeBandR2Pos','HG_FingersWholeBandR2Vel',...
    'HG_FingersSubBandsR2Pos','HG_FingersSubBandsR2Vel',...
    'HG_FingersAvgSubBandsR2Pos','HG_FingersAvgSubBandsR2Vel');

%% plot the R2 for vel and pos for delta, theta, alpha, beta

% plot the R2 for vel
figure;
set(gcf, 'Position', [300, 300, 700, 600]);
suptitle(['R2 (Velocity) for Fingers']);
plot(Delta_FingersR2Vel,'-o','color','r','LineWidth',1.5)
hold on
plot(Theta_FingersR2Vel,'-o','color','b','LineWidth',1.5)
hold on
plot(Alpha_FingersR2Vel,'-o','color','k','LineWidth',1.5)
hold on
plot(Beta_FingersR2Vel,'-o','color','m','LineWidth',1.5)
hold on
xlim([1,5])
xticks([1:5])
ylim([0,0.8])
yticks([0:0.1:0.8])
xlabel('Finger')
ylabel('R2')
set(gca,'fontsize',14)
legend('Delta','Theta','Alpha','Beta');

% plot the R2 for pos
figure;
set(gcf, 'Position', [300, 300, 700, 600]);
suptitle(['R2 (Position) for Fingers']);
plot(Delta_FingersR2Pos,'-o','color','r','LineWidth',1.5)
hold on
plot(Theta_FingersR2Pos,'-o','color','b','LineWidth',1.5)
hold on
plot(Alpha_FingersR2Pos,'-o','color','k','LineWidth',1.5)
hold on
plot(Beta_FingersR2Pos,'-o','color','m','LineWidth',1.5)
hold on
xlim([1,5])
xticks([1:5])
ylim([0,0.8])
yticks([0:0.1:0.8])
xlabel('Finger')
ylabel('R2')
set(gca,'fontsize',14)
legend('Delta','Theta','Alpha','Beta');


%% plot the R2 for vel and pos for Low Gamma

% for whole band
% plot the R2 for vel
figure;
set(gcf, 'Position', [300, 300, 700, 600]);
suptitle(['R2 (Velocity) for Fingers; WholeLG band']);
plot(LG_FingersWholeBandR2Vel.PureEnv,'-o','color','r','LineWidth',1.5)
hold on
plot(LG_FingersWholeBandR2Vel.EnvDeltaWOutMean,'-o','color','b','LineWidth',1.5)
hold on
plot(LG_FingersWholeBandR2Vel.EnvThetaWOutMean,'-o','color','k','LineWidth',1.5)
hold on
xlim([1,5])
xticks([1:5])
ylim([0,0.7])
yticks([0:0.1:0.7])
xlabel('Finger')
ylabel('R2')
set(gca,'fontsize',14)
legend('PureEnv','EnvDeltaWOutMean','EnvThetaWOutMean');

% for whole band
% plot the R2 for vel
figure;
set(gcf, 'Position', [300, 300, 700, 600]);
suptitle(['R2 (Postion) for Fingers; WholeLG band']);
plot(LG_FingersWholeBandR2Pos.PureEnv,'-o','color','r','LineWidth',1.5)
hold on
plot(LG_FingersWholeBandR2Pos.EnvDeltaWOutMean,'-o','color','b','LineWidth',1.5)
hold on
plot(LG_FingersWholeBandR2Pos.EnvThetaWOutMean,'-o','color','k','LineWidth',1.5)
hold on
xlim([1,5])
xticks([1:5])
ylim([0,0.7])
yticks([0:0.1:0.7])
xlabel('Finger')
ylabel('R2')
set(gca,'fontsize',14)

% for subband
% plot the R2 for vel
for Fi=1:5
    figure;
    set(gcf, 'Position', [300, 300, 700, 600]);
    suptitle(['R2 (Velocity) for SubBands;Finger: ',num2str(Fi)]);
    plot(LG_FingersSubBandsR2Vel(Fi).PureEnv,'-o','color','r','LineWidth',1.5)
    hold on
    plot(LG_FingersSubBandsR2Vel(Fi).EnvDeltaWOutMean,'-o','color','b','LineWidth',1.5)
    hold on
    plot(LG_FingersSubBandsR2Vel(Fi).EnvThetaWOutMean,'-o','color','k','LineWidth',1.5)
    hold on
    xlim([1,3])
    xticks([1:3])
    ylim([0,0.7])
    yticks([0:0.1:0.7])
    xlabel('Subband')
    ylabel('R2')
    set(gca,'fontsize',14)
    if Fi==1
        xlim([1,3])
        legend('PureEnv','EnvDeltaWOutMean','EnvThetaWOutMean');
    end
end

% plot the R2 for pos
for Fi=1:5
    figure;
    set(gcf, 'Position', [300, 300, 700, 600]);
    suptitle(['R2 (Position) for SubBands;Finger: ',num2str(Fi)]);
    plot(LG_FingersSubBandsR2Pos(Fi).PureEnv,'-o','color','r','LineWidth',1.5)
    hold on
    plot(LG_FingersSubBandsR2Pos(Fi).EnvDeltaWOutMean,'-o','color','b','LineWidth',1.5)
    hold on
    plot(LG_FingersSubBandsR2Pos(Fi).EnvThetaWOutMean,'-o','color','k','LineWidth',1.5)
    hold on
    xlim([1,3])
    xticks([1:3])
    ylim([0,0.7])
    yticks([0:0.1:0.7])
    xlabel('Subband')
    ylabel('R2')
    set(gca,'fontsize',14)
    if Fi==1
        xlim([1,3])
        legend('PureEnv','EnvDeltaWOutMean','EnvThetaWOutMean');
    end
end

% for the avereged bands
% plot the R2 for vel
figure;
set(gcf, 'Position', [300, 300, 700, 600]);
suptitle(['R2 (Velocity) for Fingers; Avg of SubBands']);
plot(LG_FingersAvgSubBandsR2Vel.PureEnv,'-o','color','r','LineWidth',1.5)
hold on
plot(LG_FingersAvgSubBandsR2Vel.EnvDeltaWOutMean,'-o','color','b','LineWidth',1.5)
hold on
plot(LG_FingersAvgSubBandsR2Vel.EnvThetaWOutMean,'-o','color','k','LineWidth',1.5)
xlim([1,5])
xticks([1:5])
ylim([0,0.8])
yticks([0:0.1:0.7])
xlabel('Finger')
ylabel('R2')
set(gca,'fontsize',14)
legend('PureEnv','EnvDeltaWOutMean','EnvThetaWOutMean');

% for the avereged bands
% plot the R2 for vel
figure;
set(gcf, 'Position', [300, 300, 700, 600]);
suptitle(['R2 (Position) for Fingers; Avg of SubBands']);
plot(LG_FingersAvgSubBandsR2Pos.PureEnv,'-o','color','r','LineWidth',1.5)
hold on
plot(LG_FingersAvgSubBandsR2Pos.EnvDeltaWOutMean,'-o','color','b','LineWidth',1.5)
hold on
plot(LG_FingersAvgSubBandsR2Pos.EnvThetaWOutMean,'-o','color','k','LineWidth',1.5)
xlim([1,5])
xticks([1:5])
ylim([0,0.8])
yticks([0:0.1:0.7])
xlabel('Finger')
ylabel('R2')
set(gca,'fontsize',14)
legend('PureEnv','EnvDeltaWOutMean','EnvThetaWOutMean');



%% plot the R2 for vel and pos for WholeGammaBand

% HG_FingersWholeBandR2Vel

% plot the R2 for vel
figure;
set(gcf, 'Position', [300, 300, 700, 600]);
suptitle(['R2 (Velocity) for Fingers; WholeHG band']);
plot(HG_FingersWholeBandR2Vel.EnvDeltaWOutMean,'-o','color','r','LineWidth',1.5)
hold on
plot(HG_FingersWholeBandR2Vel.EnvDeltaWithMean,'-o','color','b','LineWidth',1.5)
hold on
plot(HG_FingersWholeBandR2Vel.EnvThetaWOutMean,'-o','color','y','LineWidth',1.5)
hold on
plot(HG_FingersWholeBandR2Vel.EnvThetaWithMean,'-o','color','m','LineWidth',1.5)
hold on
plot(HG_FingersWholeBandR2Vel.EnvAlphaWOutMean,'-o','color','c','LineWidth',1.5)
hold on
plot(HG_FingersWholeBandR2Vel.EnvAlphaWithMean,'-o','color','g','LineWidth',1.5)
hold on
plot(HG_FingersWholeBandR2Vel.EnvBetaWOutMean,'-o','color','k','LineWidth',1.5)
hold on
plot(HG_FingersWholeBandR2Vel.EnvBetaWithMean,'-o','color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
xlim([1,8])
xticks([1:5])
ylim([0,0.8])
yticks([0:0.1:0.8])
xlabel('Finger')
ylabel('R2')
set(gca,'fontsize',14)
legend('EnvDeltaWOutMean','EnvDeltaWithMean','EnvThetaWOutMean','EnvThetaWithMean',...
    'EnvAlphaWOutMean','EnvAlphaWithMean','EnvBetaWOutMean','EnvBetaWithMean');

% plot the R2 for pos
figure;
set(gcf, 'Position', [300, 300, 700, 600]);
suptitle(['R2 (Postion) for Fingers; WholeHG band']);
plot(HG_FingersWholeBandR2Pos.EnvDeltaWOutMean,'-o','color','r','LineWidth',1.5)
hold on
plot(HG_FingersWholeBandR2Pos.EnvDeltaWithMean,'-o','color','b','LineWidth',1.5)
hold on
plot(HG_FingersWholeBandR2Pos.EnvThetaWOutMean,'-o','color','y','LineWidth',1.5)
hold on
plot(HG_FingersWholeBandR2Pos.EnvThetaWithMean,'-o','color','m','LineWidth',1.5)
hold on
plot(HG_FingersWholeBandR2Pos.EnvAlphaWOutMean,'-o','color','c','LineWidth',1.5)
hold on
plot(HG_FingersWholeBandR2Pos.EnvAlphaWithMean,'-o','color','g','LineWidth',1.5)
hold on
plot(HG_FingersWholeBandR2Pos.EnvBetaWOutMean,'-o','color','k','LineWidth',1.5)
hold on
plot(HG_FingersWholeBandR2Pos.EnvBetaWithMean,'-o','color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
xlim([1,8])
xticks([1:5])
ylim([0,0.8])
yticks([0:0.1:0.8])
xlabel('Finger')
ylabel('R2')
set(gca,'fontsize',14)
legend('EnvDeltaWOutMean','EnvDeltaWithMean','EnvThetaWOutMean','EnvThetaWithMean',...
    'EnvAlphaWOutMean','EnvAlphaWithMean','EnvBetaWOutMean','EnvBetaWithMean');

%% plot the R2 for vel and pos for subbands of GammaBand
% HG_FingersSubBandsR2Vel

% plot the R2 for vel
for Fi=1:5
    figure;
    set(gcf, 'Position', [300, 300, 700, 600]);
    suptitle(['R2 (Velocity) for SubBands;Finger: ',num2str(Fi)]);
    plot(HG_FingersSubBandsR2Vel(Fi).EnvDeltaWOutMean,'-o','color','r','LineWidth',1.5)
    hold on
    plot(HG_FingersSubBandsR2Vel(Fi).EnvDeltaWithMean,'-o','color','b','LineWidth',1.5)
    hold on
    plot(HG_FingersSubBandsR2Vel(Fi).EnvThetaWOutMean,'-o','color','y','LineWidth',1.5)
    hold on
    plot(HG_FingersSubBandsR2Vel(Fi).EnvThetaWithMean,'-o','color','m','LineWidth',1.5)
    hold on
    plot(HG_FingersSubBandsR2Vel(Fi).EnvAlphaWOutMean,'-o','color','c','LineWidth',1.5)
    hold on
    plot(HG_FingersSubBandsR2Vel(Fi).EnvAlphaWithMean,'-o','color','g','LineWidth',1.5)
    hold on
    plot(HG_FingersSubBandsR2Vel(Fi).EnvBetaWOutMean,'-o','color','k','LineWidth',1.5)
    hold on
    plot(HG_FingersSubBandsR2Vel(Fi).EnvBetaWithMean,'-o','color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
    xlim([1,8])
    xticks([1:8])
    ylim([0,0.8])
    yticks([0:0.1:0.8])
    xlabel('Subband')
    ylabel('R2')
    set(gca,'fontsize',14)
    if Fi==1
        xlim([1,14])
        legend('EnvDeltaWOutMean','EnvDeltaWithMean','EnvThetaWOutMean','EnvThetaWithMean',...
            'EnvAlphaWOutMean','EnvAlphaWithMean','EnvBetaWOutMean','EnvBetaWithMean');
    end
end

% plot the R2 for pos
for Fi=1:5
    figure;
    set(gcf, 'Position', [300, 300, 700, 600]);
    suptitle(['R2 (Position) for SubBands;Finger: ',num2str(Fi)]);
    plot(HG_FingersSubBandsR2Pos(Fi).EnvDeltaWOutMean,'-o','color','r','LineWidth',1.5)
    hold on
    plot(HG_FingersSubBandsR2Pos(Fi).EnvDeltaWithMean,'-o','color','b','LineWidth',1.5)
    hold on
    plot(HG_FingersSubBandsR2Pos(Fi).EnvThetaWOutMean,'-o','color','y','LineWidth',1.5)
    hold on
    plot(HG_FingersSubBandsR2Pos(Fi).EnvThetaWithMean,'-o','color','m','LineWidth',1.5)
    hold on
    plot(HG_FingersSubBandsR2Pos(Fi).EnvAlphaWOutMean,'-o','color','c','LineWidth',1.5)
    hold on
    plot(HG_FingersSubBandsR2Pos(Fi).EnvAlphaWithMean,'-o','color','g','LineWidth',1.5)
    hold on
    plot(HG_FingersSubBandsR2Pos(Fi).EnvBetaWOutMean,'-o','color','k','LineWidth',1.5)
    hold on
    plot(HG_FingersSubBandsR2Pos(Fi).EnvBetaWithMean,'-o','color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
    xlim([1,8])
    xticks([1:8])
    ylim([0,0.9])
    yticks([0:0.1:0.9])
    xlabel('Subband')
    ylabel('R2')
    set(gca,'fontsize',14)
    if Fi==1
        xlim([1,14])
        legend('EnvDeltaWOutMean','EnvDeltaWithMean','EnvThetaWOutMean','EnvThetaWithMean',...
            'EnvAlphaWOutMean','EnvAlphaWithMean','EnvBetaWOutMean','EnvBetaWithMean');
    end
end

%% plot the R2 for vel and pos for averaged subbands of GammaBand
% HG_FingersAvgSubBandsR2Vel

% plot the R2 for vel
figure;
set(gcf, 'Position', [300, 300, 700, 600]);
suptitle(['R2 (Velocity) for Fingers; Avg of SubBands']);
plot(HG_FingersAvgSubBandsR2Vel.EnvDeltaWOutMean,'-o','color','r','LineWidth',1.5)
hold on
plot(HG_FingersAvgSubBandsR2Vel.EnvDeltaWithMean,'-o','color','b','LineWidth',1.5)
hold on
plot(HG_FingersAvgSubBandsR2Vel.EnvThetaWOutMean,'-o','color','y','LineWidth',1.5)
hold on
plot(HG_FingersAvgSubBandsR2Vel.EnvThetaWithMean,'-o','color','m','LineWidth',1.5)
hold on
plot(HG_FingersAvgSubBandsR2Vel.EnvAlphaWOutMean,'-o','color','c','LineWidth',1.5)
hold on
plot(HG_FingersAvgSubBandsR2Vel.EnvAlphaWithMean,'-o','color','g','LineWidth',1.5)
hold on
plot(HG_FingersAvgSubBandsR2Vel.EnvBetaWOutMean,'-o','color','k','LineWidth',1.5)
hold on
plot(HG_FingersAvgSubBandsR2Vel.EnvBetaWithMean,'-o','color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
xlim([1,8])
xticks([1:5])
ylim([0,0.8])
yticks([0:0.1:0.8])
xlabel('Finger')
ylabel('R2')
set(gca,'fontsize',14)
legend('EnvDeltaWOutMean','EnvDeltaWithMean','EnvThetaWOutMean','EnvThetaWithMean',...
    'EnvAlphaWOutMean','EnvAlphaWithMean','EnvBetaWOutMean','EnvBetaWithMean');


% plot the R2 for pos
figure;
set(gcf, 'Position', [300, 300, 700, 600]);
suptitle(['R2 (Position) for Fingers; Avg of SubBands']);
plot(HG_FingersAvgSubBandsR2Pos.EnvDeltaWOutMean,'-o','color','r','LineWidth',1.5)
hold on
plot(HG_FingersAvgSubBandsR2Pos.EnvDeltaWithMean,'-o','color','b','LineWidth',1.5)
hold on
plot(HG_FingersAvgSubBandsR2Pos.EnvThetaWOutMean,'-o','color','y','LineWidth',1.5)
hold on
plot(HG_FingersAvgSubBandsR2Pos.EnvThetaWithMean,'-o','color','m','LineWidth',1.5)
hold on
plot(HG_FingersAvgSubBandsR2Pos.EnvAlphaWOutMean,'-o','color','c','LineWidth',1.5)
hold on
plot(HG_FingersAvgSubBandsR2Pos.EnvAlphaWithMean,'-o','color','g','LineWidth',1.5)
hold on
plot(HG_FingersAvgSubBandsR2Pos.EnvBetaWOutMean,'-o','color','k','LineWidth',1.5)
hold on
plot(HG_FingersAvgSubBandsR2Pos.EnvBetaWithMean,'-o','color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
xlim([1,8])
xticks([1:5])
ylim([0,0.8])
yticks([0:0.1:0.8])
xlabel('Finger')
ylabel('R2')
set(gca,'fontsize',14)
legend('EnvDeltaWOutMean','EnvDeltaWithMean','EnvThetaWOutMean','EnvThetaWithMean',...
    'EnvAlphaWOutMean','EnvAlphaWithMean','EnvBetaWOutMean','EnvBetaWithMean');

%% plot the R2 for vel and pos for all comparison

% plot the R2 for vel
for Fi=1:5
    figure;
    set(gcf, 'Position', [300, 300, 700, 600]);
    suptitle(['R2 (Velocity) for Finger: ',num2str(Fi)]);
    value1=HG_FingersAvgSubBandsR2Vel.EnvDeltaWOutMean(Fi);
    plot(0:1:15,ones(1,16)*value1,'r','LineWidth',1.5)
    hold on
    value2=HG_FingersWholeBandR2Vel.EnvDeltaWOutMean(Fi);
    plot(0:1:15,ones(1,16)*value2,'b','LineWidth',1.5)
    hold on
    plot(HG_FingersSubBandsR2Vel(Fi).EnvDeltaWOutMean,'-o','color','k','LineWidth',1.5)
    hold on
    plot(HG_FingersSubBandsR2Vel(Fi).EnvThetaWOutMean,'-o','color','m','LineWidth',1.5)
    hold on
    xlim([1,8])
    xticks([1:8])
    ylim([0.4,0.8])
    yticks([0.4:0.1:0.8])
    xlabel('Subband')
    ylabel('R2')
    set(gca,'fontsize',14)
    
    if Fi==1
        xlim([1,8])
        legend('AvgSubBandsDelta','WholeHGBandDelta','SubBandsDelta','SubBandsTheta');
    end
end

% plot the R2 for pos
for Fi=1:5
    figure;
    set(gcf, 'Position', [300, 300, 700, 600]);
    suptitle(['R2 (Position) for Finger: ',num2str(Fi)]);
    value1=HG_FingersAvgSubBandsR2Pos.EnvDeltaWOutMean(Fi);
    plot(0:1:15,ones(1,16)*value1,'r','LineWidth',1.5)
    hold on
    value2=HG_FingersWholeBandR2Pos.EnvDeltaWOutMean(Fi);
    plot(0:1:15,ones(1,16)*value2,'b','LineWidth',1.5)
    hold on
    plot(HG_FingersSubBandsR2Pos(Fi).EnvDeltaWOutMean,'-o','color','k','LineWidth',1.5)
    hold on
    plot(HG_FingersSubBandsR2Pos(Fi).EnvThetaWOutMean,'-o','color','m','LineWidth',1.5)
    hold on
    xlim([1,8])
    xticks([1:8])
    ylim([0.4,0.9])
    yticks([0.4:0.1:0.9])
    xlabel('Subband')
    ylabel('R2')
    set(gca,'fontsize',14)
    
    if Fi==1
        xlim([1,8])
        legend('AvgSubBandsDelta','WholeHGBandDelta','SubBandsDelta','SubBandsTheta');
    end
end

%% plot the R2 for vel and pos; Comparison: Pure Delta & Delta of Avg-HG

% plot the R2 for vel
figure;
set(gcf, 'Position', [300, 300, 700, 600]);
suptitle(['R2 (Velocity) for Fingers']);
plot(HG_FingersAvgSubBandsR2Vel.EnvDeltaWOutMean,'-o','color','r','LineWidth',1.5)
hold on
plot(Delta_FingersR2Vel,'-o','color','b','LineWidth',1.5)
hold on

xlim([1,5])
xticks([1:5])
ylim([0,0.8])
yticks([0:0.1:0.8])
xlabel('Finger')
ylabel('R2')
set(gca,'fontsize',14)
legend('Delta of Avg-HG','Pure Delta');

% plot the R2 for pos
figure;
set(gcf, 'Position', [300, 300, 700, 600]);
suptitle(['R2 (Position) for Fingers']);
plot(HG_FingersAvgSubBandsR2Pos.EnvDeltaWOutMean,'-o','color','r','LineWidth',1.5)
hold on
plot(Delta_FingersR2Pos,'-o','color','b','LineWidth',1.5)
hold on

xlim([1,5])
xticks([1:5])
ylim([0,0.8])
yticks([0:0.1:0.8])
xlabel('Finger')
ylabel('R2')
set(gca,'fontsize',14)
legend('Delta of Avg-HG','Pure Delta');

%% Calculating the weights of MLR fitting for delta and HG
Kin_Pos={Posq_F1(:,1) Posq_F2(:,3) Posq_F3(:,3) Posq_F4(:,3) Posq_F5(:,3)};
Kin_Vel={AbsVel_F1 AbsVel_F2 AbsVel_F3 AbsVel_F4 AbsVel_F5};

% Calculating the weights for pure delta
% prediction the velocities and positions and cacculating the R2

for Fi=1:5
    Figures=0;
    WeightVel=[];
    WeightPos=[];
    Input_dynamics=Delta_Fingers(Fi).Hilbert;
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    WeightVel=MultipleRegFuncCalWeight(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    Delta_WeightsVel(Fi).Finger=WeightVel;
    WeightPos=MultipleRegFuncCalWeight(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    Delta_WeightsPos(Fi).Finger=WeightPos;

end

% save weights for delta
save('Delta_AllWeights.mat','Delta_WeightsVel','Delta_WeightsPos');

% Calculating the weights for avereged high gamma
for Fi=1:5
    Figures=0;
    WeightVel=[];
    WeightPos=[];
    Input_dynamics=abs(hilbert(HG_FingersAvgSubBands(Fi).EnvDeltaWOutMean));
    Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
    Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
    WeightVel=MultipleRegFuncCalWeight(Fi,Input_dynamics,Output_DynamicsVel,Figures);
    HG_WeightsVel(Fi).Finger=WeightVel;
    WeightPos=MultipleRegFuncCalWeight(Fi,Input_dynamics,Output_DynamicsPos,Figures);
    HG_WeightsPos(Fi).Finger=WeightPos;
end

% save weights for delta
save('HG_AllWeights.mat','HG_WeightsVel','HG_WeightsPos');


%% Ignore the following for right now: Compute the R2 for indiviudal Channel and observe the grid map

Kin_Pos={Posq_F1(:,1) Posq_F2(:,3) Posq_F3(:,3) Posq_F4(:,3) Posq_F5(:,3)};
Kin_Vel={AbsVel_F1 AbsVel_F2 AbsVel_F3 AbsVel_F4 AbsVel_F5};

% for HG_DeltaofAvgHG 
R2ChVel=[];
R2ChPos=[];
for Fi=5
    Figures=0;
    Input_dynamics=abs(hilbert(HG_DeltaofAvgHG(Fi).Finger));
    for ch=1:256
        Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
        Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
        R2Vel=MultipleRegFunc(Fi,Input_dynamics(:,ch),Output_DynamicsVel,Figures);
        R2ChVel=[R2ChVel;R2Vel];
        R2Pos=MultipleRegFunc(Fi,Input_dynamics(:,ch),Output_DynamicsPos,Figures);
        R2ChPos=[R2ChPos;R2Pos];
    end
end

Weights=R2ChVel;
Nch_Record=256;
% grid layout
Ch_num_1=1:Nch_Record;
Ch_num_2=reshape(Ch_num_1,[16,16]);
% Ch_num_3 is the grid layout
Ch_num_3=rot90(rot90(Ch_num_2));

[R,C]=size(Ch_num_3);

% for vel for high gamma 
for Fi=5
    Layout=zeros(R,C);
    for r=1:R
        for c=1:C
            Layout(r,c)=Weights(Ch_num_3(r,c),1);
        end
    end
    figure;
    imagesc(Layout)
    title(['Finger: ',num2str(Fi)])
    colorbar
    
end

%% Make a movie for on a grid showing the change of Hilberts power for delta and HG































