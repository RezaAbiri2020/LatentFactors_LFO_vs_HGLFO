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

%% loading and breaking raw ECoG data into trials
load('E:\ECoGLeapMotion\DataPatientOne\ECoGData\EC169_ecog_kinematics_25Hz.mat');

if 0
    for i=193:256
        subplot(8,8,i-192)
        plot(raw_fullFs{1,1}(1:1000,i))
        
    end
    
    figure;
    plot(raw_fullFs{1,1}(1:1000,64+49))
        
end

% Final list of bad channels by observation
BadChs=[64+1, 64+2, 64+28, 64+49, 64+62, 64+64,...
    128+1, 128+2, 128+16, 128+20, 128+24, 128+25, 128+26, 128+30, 128+31, 128+32, 128+34, 128+42, 128+46, 128+64,...
    192+1];

All_Index=ones(size(raw_fullFs{1,1},2),1);
All_Index(BadChs')=0;
Selected_Chs=logical(All_Index(1:256,1));

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

%% 1- generating the max pos and max vel using the dominant dimension

%% 1-1- F1 movement analysis: Finding the peak points for transition from flextion to extention or extension to flextion

%F1 check up
j=1;
for i=1:5
    figure
    plot(fingers{1,1}(j).joints(i).spos)
end
% the position plot for end point 
figure;
plot(fingers{1,1}(1).joints(5).spos)
X1 = fingers{1,1}(1).joints(5).spos(:,1);
X2 = fingers{1,1}(1).joints(5).spos(:,2);
X3 = fingers{1,1}(1).joints(5).spos(:,3);
figure;
plot(((X1-X1(1)).^2+(X2-X2(1)).^2+(X3-X3(1)).^2).^0.5)

Xq_F1=0:1/(Fs):(time{1,1}(end)-time{1,1}(1));
Posq_F1=interp1(time{1,1}-(time{1,1}(1)*ones(size(time{1,1},1),1)),fingers{1,1}(1).joints(5).spos,Xq_F1);
X1 = Posq_F1(:,1);
X2 = Posq_F1(:,2);
X3 = Posq_F1(:,3);
figure;
plot(((X1-X1(1)).^2+(X2-X2(1)).^2+(X3-X3(1)).^2).^0.5)

half_window=500; %points; about 1 sec: for finding local max and local min in position data

% for finger 1 
[F1_n,F1_m]=size(Xq_F1');

% max and min sample points
F1_E2F_sample=[];
F1_F2E_sample=[];

% finding max and min time stamps
for i=15*half_window:F1_n-(half_window+1)
    
    logic1=(Posq_F1(i-half_window:i+half_window,1)< Posq_F1(i,1));
    
    if sum(logic1)==2*half_window
        F1_E2F_sample=[F1_E2F_sample;i];
        
    end
    
    logic2=(Posq_F1(i-half_window:i+half_window,1)> Posq_F1(i,1));
    
    if sum(logic2)==2*half_window
        F1_F2E_sample=[F1_F2E_sample;i];
        
    end
    
    
end

F1_E2F_sample=[3221;6750;F1_E2F_sample];
F1_F2E_sample=[4685;F1_F2E_sample([1:8])];

% check on the plot
figure;
plot(Posq_F1(:,1))
hold on; plot(F1_F2E_sample,30,'og')
hold on; plot(F1_E2F_sample,-15,'or')


%% 1-2- F1 movement analysis: Finding the max velocity during transition from flextion to extention or extension to flextion
% derivative of original signal and not Posq_F1

Vel_F1=diff(fingers{1,1}(1).joints(5).spos);
Velq_F1=interp1(time{1,1}-(time{1,1}(1)*ones(size(time{1,1},1),1)),[Vel_F1;0,0,0],Xq_F1);
DimVel_F1=Velq_F1(:,1);

Index_F1=sort([1; F1_E2F_sample; F1_F2E_sample; length(Posq_F1)]);

F1_E2F_Vels=[];

for i=1:2:(length(Index_F1)-1)
    Section=DimVel_F1(Index_F1(i):(Index_F1(i+1)));
    [m_S n_S]=max(Section);
    F1_E2F_Vels=[F1_E2F_Vels; (n_S+Index_F1(i)-1)];  
    
end

F1_F2E_Vels=[];

for i=2:2:(length(Index_F1)-1)
    Section=DimVel_F1(Index_F1(i):(Index_F1(i+1)));
    [m_S n_S]=min(Section);
    F1_F2E_Vels=[F1_F2E_Vels; (n_S+Index_F1(i)-1)];
    
end

% check on the plot
figure;
plot(DimVel_F1)
hold on; plot(F1_F2E_Vels,-2,'ok')
hold on; plot(F1_E2F_Vels,2,'ob')

plot(Posq_F1(:,1))
hold on; plot(F1_F2E_sample,30,'og')
hold on; plot(F1_E2F_sample,-15,'or')
hold on; vline(F1_F2E_Vels,'g--')
hold on; vline(F1_E2F_Vels,'r--')


FingersKinInfo.Finger(1).Fs_MaxPos=F1_E2F_sample;
FingersKinInfo.Finger(1).Es_MaxPos=F1_F2E_sample;

FingersKinInfo.Finger(1).Fs_MaxVel = F1_E2F_Vels; 
FingersKinInfo.Finger(1).Es_MaxVel = F1_F2E_Vels;

FingersKinData.Finger(1).PurePos=Posq_F1(:,1);
FingersKinData.Finger(1).PureVel=DimVel_F1;


%% 1-3 F2 movement analysis: Finding the peak points for transition from flextion to extention or extension to flextion
% F2 check up

j=2;
for i=1:5
    figure
    plot(fingers{1,2}(j).joints(i).spos)
end

% the position plot for end point 
figure;
plot(fingers{1,2}(2).joints(5).spos)

Xq_F2=0:1/(Fs):(time{1,2}(end)-time{1,2}(1));
Posq_F2=interp1(time{1,2}-(time{1,2}(1)*ones(size(time{1,2},1),1)),fingers{1,2}(2).joints(5).spos,Xq_F2);
%plot(Posq_F2)
X1 = Posq_F2(:,1);
X2 = Posq_F2(:,2);
X3 = Posq_F2(:,3);
figure;
plot(((X1-X1(1)).^2+(X2-X2(1)).^2+(X3-X3(1)).^2).^0.5)


half_window=300; %points; about 1 sec: for finding local max and local min in position data

% for finger 1 
[F2_n,F2_m]=size(Xq_F2');

% max and min sample points
F2_F2E_sample=[];
F2_E2F_sample=[];

% finding max and min time stamps
for i=1.5*half_window:F2_n-(half_window+1)
    
    logic1=(Posq_F2(i-half_window:i+half_window,3)< Posq_F2(i,3));
    
    if sum(logic1)==2*half_window
        F2_E2F_sample=[F2_E2F_sample;i];
        
    end
    
    logic2=(Posq_F2(i-half_window:i+half_window,3)> Posq_F2(i,3));
    
    if sum(logic2)==2*half_window
        F2_F2E_sample=[F2_F2E_sample;i];
        
    end
     
end

% check on the plot
figure;
plot(Posq_F2(:,3))
hold on; plot(F2_F2E_sample,30,'og')
hold on; plot(F2_E2F_sample,0,'or')

%% 1-4 F2 movement analysis: Finding the max velocity during transition from flextion to extention or extension to flextion

Vel_F2=diff(fingers{1,2}(2).joints(5).spos);
Velq_F2=interp1(time{1,2}-(time{1,2}(1)*ones(size(time{1,2},1),1)),[Vel_F2;0,0,0],Xq_F2);
DimVel_F2=Velq_F2(:,3);

Index_F2=sort([1; F2_E2F_sample; F2_F2E_sample]);

F2_E2F_Vels=[];

for i=1:2:(length(Index_F2)-1)
    Section=DimVel_F2(Index_F2(i):(Index_F2(i+1)));
    [m_S n_S]=max(Section);
    F2_E2F_Vels=[F2_E2F_Vels; (n_S+Index_F2(i)-1)];  
    
end

F2_F2E_Vels=[];

for i=2:2:(length(Index_F2)-1)
    Section=DimVel_F2(Index_F2(i):(Index_F2(i+1)));
    [m_S n_S]=min(Section);
    F2_F2E_Vels=[F2_F2E_Vels; (n_S+Index_F2(i)-1)];
    
end

% check on the plot
figure;
plot(DimVel_F2)
hold on; plot(F2_F2E_Vels,6,'ok')
hold on; plot(F2_E2F_Vels,10,'ob')

figure;
hold on;
plot(Posq_F2(:,3))
hold on; plot(F2_F2E_sample,30,'og')
hold on; plot(F2_E2F_sample,0,'or')
hold on; vline(F2_F2E_Vels,'g--')
hold on; vline(F2_E2F_Vels,'r--')

FingersKinInfo.Finger(2).Fs_MaxPos=F2_E2F_sample;
FingersKinInfo.Finger(2).Es_MaxPos=F2_F2E_sample;

FingersKinInfo.Finger(2).Fs_MaxVel = F2_E2F_Vels; 
FingersKinInfo.Finger(2).Es_MaxVel = F2_F2E_Vels;

FingersKinData.Finger(2).PurePos=Posq_F2(:,3);
FingersKinData.Finger(2).PureVel=DimVel_F2;


%% 1-5 F3 movement analysis: Finding the peak points for transition from flextion to extention or extension to flextion

% F3 check up
j=3;
for i=1:5
    figure
    plot(fingers{1,3}(j).joints(i).spos)
end

% the position plot for end point 
figure;
plot(fingers{1,3}(3).joints(5).spos)

Xq_F3=0:1/(Fs):(time{1,3}(end)-time{1,3}(1));
Posq_F3=interp1(time{1,3}-(time{1,3}(1)*ones(size(time{1,3},1),1)),fingers{1,3}(3).joints(5).spos,Xq_F3);
%plot(Posq_F2)

half_window=350; %points; about 1 sec: for finding local max and local min in position data

% for finger 1 
[F3_n,F3_m]=size(Xq_F3');

% max and min sample points
F3_F2E_sample=[];
F3_E2F_sample=[];

% finding max and min time stamps
for i=5*half_window:F3_n-(half_window+1)
    
    logic1=(Posq_F3(i-half_window:i+half_window,3)< Posq_F3(i,3));
    
    if sum(logic1)==2*half_window
        F3_E2F_sample=[F3_E2F_sample;i];
        
    end
    
    logic2=(Posq_F3(i-half_window:i+half_window,3)> Posq_F3(i,3));
    
    if sum(logic2)==2*half_window
        F3_F2E_sample=[F3_F2E_sample;i];
        
    end
    
    
end

%F3_F2E_sample=F3_F2E_sample([1:3,5:end]);
% check on the plot
figure;
plot(Posq_F3(:,3))
hold on; plot(F3_F2E_sample,20,'og')
hold on; plot(F3_E2F_sample,-20,'or')


%% 1-6 F3 movement analysis: Finding the max velocity during transition from flextion to extention or extension to flextion

Vel_F3=diff(fingers{1,3}(3).joints(5).spos);
Velq_F3=interp1(time{1,3}-(time{1,3}(1)*ones(size(time{1,3},1),1)),[Vel_F3;0,0,0],Xq_F3);
DimVel_F3=Velq_F3(:,3);
%plot(AbsVel_F1)

Index_F3=sort([1; F3_E2F_sample; F3_F2E_sample]);

F3_E2F_Vels=[];

for i=1:2:(length(Index_F3)-1)
    Section=DimVel_F3(Index_F3(i):(Index_F3(i+1)));
    [m_S n_S]=max(Section);
    F3_E2F_Vels=[F3_E2F_Vels; (n_S+Index_F3(i)-1)];  
    
end

F3_F2E_Vels=[];

for i=2:2:(length(Index_F3)-1)
    Section=DimVel_F3(Index_F3(i):(Index_F3(i+1)));
    [m_S n_S]=min(Section);
    F3_F2E_Vels=[F3_F2E_Vels; (n_S+Index_F3(i)-1)];
    
end


% check on the plot
figure;
plot(DimVel_F3)
hold on; plot(F3_F2E_Vels,6,'ok')
hold on; plot(F3_E2F_Vels,10,'ob')

figure;
plot(Posq_F3(:,3))
hold on; plot(F3_F2E_sample,20,'og')
hold on; plot(F3_E2F_sample,-20,'or')
hold on; vline(F3_F2E_Vels,'g--')
hold on; vline(F3_E2F_Vels,'r--')

FingersKinInfo.Finger(3).Fs_MaxPos=F3_E2F_sample;
FingersKinInfo.Finger(3).Es_MaxPos=F3_F2E_sample;

FingersKinInfo.Finger(3).Fs_MaxVel = F3_E2F_Vels; 
FingersKinInfo.Finger(3).Es_MaxVel = F3_F2E_Vels;

FingersKinData.Finger(3).PurePos=Posq_F3(:,3);
FingersKinData.Finger(3).PureVel=DimVel_F3;

%% 1-7 F4 movement analysis: Finding the peak points for transition from flextion to extention or extension to flextion

% F4 check up
j=4;
for i=1:5
    figure
    plot(fingers{1,4}(j).joints(i).spos)
end

% the position plot for end point 
figure;
plot(fingers{1,4}(4).joints(5).spos)

Xq_F4=0:1/(Fs):(time{1,4}(end)-time{1,4}(1));
Posq_F4=interp1(time{1,4}-(time{1,4}(1)*ones(size(time{1,4},1),1)),fingers{1,4}(4).joints(4).spos,Xq_F4);
plot(Posq_F4(:,3))

F4_E2F_sample=[18900;23080;26000;27650;29330;31050;34950;36830;38740];
F4_F2E_sample=[20020;24490;26930;28640;30310;32110;36100;37930;39460];


% check on the plot
figure;
plot(Posq_F4(:,3))
hold on; plot(F4_F2E_sample,10,'og')
hold on; plot(F4_E2F_sample,-10,'or')

%% 1-8 F4 movement analysis: Finding the max velocity during transition from flextion to extention or extension to flextion

Vel_F4=diff(fingers{1,4}(4).joints(5).spos);
Velq_F4=interp1(time{1,4}-(time{1,4}(1)*ones(size(time{1,4},1),1)),[Vel_F4;0,0,0],Xq_F4);
DimVel_F4=Velq_F4(:,3);

Index_F4=sort([17610; F4_E2F_sample; F4_F2E_sample]);

F4_E2F_Vels=[];

for i=1:2:(length(Index_F4)-1)
    Section=DimVel_F4(Index_F4(i):(Index_F4(i+1)));
    [m_S n_S]=max(Section);
    F4_E2F_Vels=[F4_E2F_Vels; (n_S+Index_F4(i)-1)];  
    
end

F4_F2E_Vels=[];

for i=2:2:(length(Index_F4)-1)
    Section=DimVel_F4(Index_F4(i):(Index_F4(i+1)));
    [m_S n_S]=min(Section);
    F4_F2E_Vels=[F4_F2E_Vels; (n_S+Index_F4(i)-1)];
    
end

F4_E2F_Vels=[F4_E2F_Vels(1);22420;F4_E2F_Vels(3:end)];


% check on the plot
figure;
plot(DimVel_F4)
hold on; plot(F4_F2E_Vels,6,'ok')
hold on; plot(F4_E2F_Vels,8,'ob')

figure;
plot(Posq_F4(:,3))
hold on; plot(F4_F2E_sample,10,'og')
hold on; plot(F4_E2F_sample,-10,'or')
hold on; vline(F4_F2E_Vels,'g--')
hold on; vline(F4_E2F_Vels,'r--')

FingersKinInfo.Finger(4).Fs_MaxPos=F4_E2F_sample;
FingersKinInfo.Finger(4).Es_MaxPos=F4_F2E_sample;

FingersKinInfo.Finger(4).Fs_MaxVel = F4_E2F_Vels; 
FingersKinInfo.Finger(4).Es_MaxVel = F4_F2E_Vels;

FingersKinData.Finger(4).PurePos=Posq_F4(:,3);
FingersKinData.Finger(4).PureVel=DimVel_F4;


%% 1-9 F5 movement analysis: Finding the peak points for transition from flextion to extention or extension to flextion

% F4 check up
j=5;
for i=1:5
    figure
    plot(fingers{1,5}(j).joints(i).spos)
end
% the position plot for end point 
figure;
plot(fingers{1,5}(5).joints(5).spos)

Xq_F5=0:1/(Fs):(time{1,5}(end)-time{1,5}(1));
Posq_F5=interp1(time{1,5}-(time{1,5}(1)*ones(size(time{1,5},1),1)),fingers{1,5}(5).joints(5).spos,Xq_F5);
plot(Posq_F5(:,3))


half_window=350; %points; about 1 sec: for finding local max and local min in position data

% for finger 1 
[F5_n,F5_m]=size(Xq_F5');

% max and min sample points
F5_F2E_sample=[];
F5_E2F_sample=[];

% finding max and min time stamps
for i=18160:F5_n-(half_window+1)
    
    logic1=(Posq_F5(i-half_window:i+half_window,3)< Posq_F5(i,3));
    
    if sum(logic1)==2*half_window
        F5_E2F_sample=[F5_E2F_sample;i];
        
    end
    
    logic2=(Posq_F5(i-half_window:i+half_window,3)> Posq_F5(i,3));
    
    if sum(logic2)==2*half_window
        F5_F2E_sample=[F5_F2E_sample;i];
        
    end
    
    
end

% check on the plot
figure;
plot(Posq_F5(:,3))
hold on; plot(F5_F2E_sample,10,'og')
hold on; plot(F5_E2F_sample,-20,'or')

%% 1-10 F5 movement analysis: Finding the max velocity during transition from flextion to extention or extension to flextion

Vel_F5=diff(fingers{1,5}(5).joints(5).spos);
Velq_F5=interp1(time{1,5}-(time{1,5}(1)*ones(size(time{1,5},1),1)),[Vel_F5;0,0,0],Xq_F5);
DimVel_F5=Velq_F5(:,3);

Index_F5=sort([18160; F5_E2F_sample; F5_F2E_sample]);

F5_E2F_Vels=[];

for i=1:2:(length(Index_F5)-1)
    Section=DimVel_F5(Index_F5(i):(Index_F5(i+1)));
    [m_S n_S]=max(Section);
    F5_E2F_Vels=[F5_E2F_Vels; (n_S+Index_F5(i)-1)];  
    
end

F5_F2E_Vels=[];

for i=2:2:(length(Index_F5)-1)
    Section=DimVel_F5(Index_F5(i):(Index_F5(i+1)));
    [m_S n_S]=min(Section);
    F5_F2E_Vels=[F5_F2E_Vels; (n_S+Index_F5(i)-1)];
    
end


% check on the plot
figure;
plot(DimVel_F5)
hold on; plot(F5_F2E_Vels,6,'ok')
hold on; plot(F5_E2F_Vels,10,'ob')

figure;
plot(Posq_F5(:,3))
hold on; plot(F5_F2E_sample,10,'og')
hold on; plot(F5_E2F_sample,-20,'or')
hold on; vline(F5_F2E_Vels,'g--')
hold on; vline(F5_E2F_Vels,'r--')

FingersKinInfo.Finger(5).Fs_MaxPos=F5_E2F_sample;
FingersKinInfo.Finger(5).Es_MaxPos=F5_F2E_sample;

FingersKinInfo.Finger(5).Fs_MaxVel = F5_E2F_Vels; 
FingersKinInfo.Finger(5).Es_MaxVel = F5_F2E_Vels;

FingersKinData.Finger(5).PurePos=Posq_F5(:,3);
FingersKinData.Finger(5).PureVel=DimVel_F5;

%% Referencing options for all Fingers
Modification=[length(Xq_F1),length(Xq_F2),length(Xq_F3),length(Xq_F4),length(Xq_F5)];
Nch_Record = 256;
for Fi=1:5
    % z score modification to sync
    ECoGdata = zscore(raw_fullFs{1,Fi});
    ECoGdata_Modified = ECoGdata(1:Modification(Fi),:);
    ECoGdata_1 = ECoGdata_Modified';

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
save('E:\ECoGLeapMotion\ResultsGroupAnalysis\github_Branch_V3/Subject1.mat',...
    'ECoG_data','FingersKinData','FingersKinInfo','Fs');

