

% Extract trials around a window of ~1s for max vel for flexion and for extension

close all;
clear all;
clc;

%% loading the data 

% subject 2: load for ECoG 25Hz
load('/media/reza/WindowsDrive/ECoGLeapMotion/DataPatientTwo/AllData/ecog_kinematics_25Hz.mat')
clear fingers hand

% subject 2: load for kinematics 25Hz
load('/media/reza/WindowsDrive/ECoGLeapMotion/DataPatientTwo/AllData/kinData_EC171_25Hz.mat')


Fs = 25;


%% calculate the PC for kinematics per finger and per grasp 

for Movement = 1:8 
    
    Kin_Data = [];
    
    for Fi = 1:5
        
        for joint = 1:5
            
            Kin_Data = [ Kin_Data, fingers{1, Movement}(Fi).joints(joint).spos];
        end 
            
    end
    
    % PCA per movement
    [coeff,score,latent,tsquared,explained] = pca(Kin_Data);
    %plot(score(:,1),'b')
    %hold on
    %plot(diff(score(:,1)),'r')
    
    PC_Pos(Movement).Movement = score(:,1);
    Vel = [0; diff(score(:,1))];
    PC_Vel(Movement).Movement = Vel;   
    
end 

%% the bumps inside PCs should be in sequence of flexion and extension
% becasue of difference in referece the sign of PC is changing
% making all similar in bump sequencing 

PC_Pos(1).Movement = PC_Pos(1).Movement;
PC_Pos(2).Movement = PC_Pos(2).Movement;
PC_Pos(3).Movement = PC_Pos(3).Movement;
PC_Pos(4).Movement = PC_Pos(4).Movement;
PC_Pos(5).Movement = -PC_Pos(5).Movement;
PC_Pos(6).Movement = -PC_Pos(6).Movement;
PC_Pos(7).Movement = PC_Pos(7).Movement;
PC_Pos(8).Movement = PC_Pos(8).Movement;

PC_Vel(1).Movement = PC_Vel(1).Movement;
PC_Vel(2).Movement = PC_Vel(2).Movement;
PC_Vel(3).Movement = PC_Vel(3).Movement;
PC_Vel(4).Movement = PC_Vel(4).Movement;
PC_Vel(5).Movement = -PC_Vel(5).Movement;
PC_Vel(6).Movement = -PC_Vel(6).Movement;
PC_Vel(7).Movement = PC_Vel(7).Movement;
PC_Vel(8).Movement = PC_Vel(8).Movement;


%% calculate the window that flexion and extension is happening using position data
 % then calcualte indices for max vel for flexion and extension 

%% Finger 1 position
half_window=37; %points; about 1 sec: for finding local max and local min in position data

% for finger 1 
Movement = 1;
Posq_F1 = PC_Pos(Movement).Movement;
[F1_n,F1_m]=size(Posq_F1);

% finding max and min time stamps
j=1;
jj=1;
for i=2*half_window:F1_n-(3*half_window+1)
    
    logic1=(Posq_F1(i-half_window:i+half_window,1)> Posq_F1(i,1));
    
    if sum(logic1)==2*half_window
        FingersKinInfo.Finger(1).Fs_MaxPos(j,1)=i;
        j=j+1;
    end
    
    logic2=(Posq_F1(i-half_window:i+half_window,1)< Posq_F1(i,1));
    
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


%% Finger 1 velocity

Index_F1=sort([1; FingersKinInfo.Finger(1).Fs_MaxPos;...
    FingersKinInfo.Finger(1).Es_MaxPos; F1_n]);

AbsVel_F1 = abs(PC_Vel(1).Movement); 

j=1;
for i=1:2:(length(Index_F1)-1)
    Section=AbsVel_F1(Index_F1(i):(Index_F1(i+1)));
    [m_S n_S]=max(Section);
    FingersKinInfo.Finger(1).Fs_MaxVel(j,1)=(n_S+Index_F1(i)-1);  
    j=j+1;
end

j=1;
for i=2:2:(length(Index_F1)-1)
    Section=AbsVel_F1(Index_F1(i):(Index_F1(i+1)));
    [m_S n_S]=max(Section);
    FingersKinInfo.Finger(1).Es_MaxVel(j,1)=(n_S+Index_F1(i)-1);
    j=j+1;
end


% check on the plot
figure;
plot(AbsVel_F1)
hold on; plot(FingersKinInfo.Finger(1).Fs_MaxVel,10,'ok')
hold on; plot(FingersKinInfo.Finger(1).Es_MaxVel,10,'ob')

figure;
plot(Posq_F1(:,1))
hold on; plot(FingersKinInfo.Finger(1).Fs_MaxPos,-30,'og')
hold on; plot(FingersKinInfo.Finger(1).Es_MaxPos,-5,'or')
hold on; vline(FingersKinInfo.Finger(1).Fs_MaxVel,'g--')
hold on; vline(FingersKinInfo.Finger(1).Es_MaxVel,'r--')

FingersKinData.Finger(1).PurePos=Posq_F1(:,1);
FingersKinData.Finger(1).PureVel=PC_Vel(1).Movement;


%% Finger 2 position
half_window=37; %points; about 1 sec: for finding local max and local min in position data

% for finger 2 
Movement = 2;
Posq_F2 = PC_Pos(Movement).Movement;
[F2_n,F2_m]=size(Posq_F2);

% finding max and min time stamps
j=1;
jj=1;
for i=2*half_window+10:F2_n-(half_window+1)
    
    logic1=(Posq_F2(i-half_window:i+half_window,1)> Posq_F2(i,1));
    
    if sum(logic1)==2*half_window
        FingersKinInfo.Finger(2).Fs_MaxPos(j,1)=i;
        j=j+1;
    end
    
    logic2=(Posq_F2(i-half_window:i+half_window,1)< Posq_F2(i,1));
    
    if sum(logic2)==2*half_window
        FingersKinInfo.Finger(2).Es_MaxPos(jj,1)=i;
        jj=jj+1;
    end
    
    
end

% check on the plot
figure;
plot(Posq_F2(:,1))
hold on; plot(FingersKinInfo.Finger(2).Fs_MaxPos,-10,'og')
hold on; plot(FingersKinInfo.Finger(2).Es_MaxPos,-5,'or')

%% Finger 2 Velocity

Index_F2=sort([1; FingersKinInfo.Finger(2).Fs_MaxPos;...
    FingersKinInfo.Finger(2).Es_MaxPos; F2_n]);

AbsVel_F2 = abs(PC_Vel(2).Movement); 

j=1;
for i=1:2:(length(Index_F2)-1)
    Section=AbsVel_F2(Index_F2(i):(Index_F2(i+1)));
    [m_S n_S]=max(Section);
    FingersKinInfo.Finger(2).Fs_MaxVel(j,1)=(n_S+Index_F2(i)-1);  
    j=j+1;
end

j=1;
for i=2:2:(length(Index_F2)-1)
    Section=AbsVel_F2(Index_F2(i):(Index_F2(i+1)));
    [m_S n_S]=max(Section);
    FingersKinInfo.Finger(2).Es_MaxVel(j,1)=(n_S+Index_F2(i)-1);
    j=j+1;
end


% check on the plot
figure;
plot(AbsVel_F2)
hold on; plot(FingersKinInfo.Finger(2).Fs_MaxVel,10,'ok')
hold on; plot(FingersKinInfo.Finger(2).Es_MaxVel,10,'ob')

figure;
plot(Posq_F2(:,1))
hold on; plot(FingersKinInfo.Finger(2).Fs_MaxPos,-30,'og')
hold on; plot(FingersKinInfo.Finger(2).Es_MaxPos,-5,'or')
hold on; vline(FingersKinInfo.Finger(2).Fs_MaxVel,'g--')
hold on; vline(FingersKinInfo.Finger(2).Es_MaxVel,'r--')

FingersKinData.Finger(2).PurePos=Posq_F2(:,1);
FingersKinData.Finger(2).PureVel=PC_Vel(2).Movement;

%% Finger 3 position
half_window=37; %points; about 1 sec: for finding local max and local min in position data

% for finger 3 
Movement = 3;
Posq_F3 = PC_Pos(Movement).Movement;
[F3_n,F3_m]=size(Posq_F3);

% finding max and min time stamps
j=1;
jj=1;
for i=2*half_window+10:F3_n-(half_window+1)
    
    logic1=(Posq_F3(i-half_window:i+half_window,1)> Posq_F3(i,1));
    
    if sum(logic1)==2*half_window
        FingersKinInfo.Finger(3).Fs_MaxPos(j,1)=i;
        j=j+1;
    end
    
    logic2=(Posq_F3(i-half_window:i+half_window,1)< Posq_F3(i,1));
    
    if sum(logic2)==2*half_window
        FingersKinInfo.Finger(3).Es_MaxPos(jj,1)=i;
        jj=jj+1;
    end
    
    
end

% check on the plot
figure;
plot(Posq_F3(:,1))
hold on; plot(FingersKinInfo.Finger(3).Fs_MaxPos,-10,'og')
hold on; plot(FingersKinInfo.Finger(3).Es_MaxPos,-5,'or')

%% Finger 3 Velocity

Index_F3=sort([1; FingersKinInfo.Finger(3).Fs_MaxPos;...
    FingersKinInfo.Finger(3).Es_MaxPos; F3_n]);

AbsVel_F3 = abs(PC_Vel(3).Movement); 

j=1;
for i=1:2:(length(Index_F3)-1)
    Section=AbsVel_F3(Index_F3(i):(Index_F3(i+1)));
    [m_S n_S]=max(Section);
    FingersKinInfo.Finger(3).Fs_MaxVel(j,1)=(n_S+Index_F3(i)-1);  
    j=j+1;
end

j=1;
for i=2:2:(length(Index_F3)-1)
    Section=AbsVel_F3(Index_F3(i):(Index_F3(i+1)));
    [m_S n_S]=max(Section);
    FingersKinInfo.Finger(3).Es_MaxVel(j,1)=(n_S+Index_F3(i)-1);
    j=j+1;
end


% check on the plot
figure;
plot(AbsVel_F3)
hold on; plot(FingersKinInfo.Finger(3).Fs_MaxVel,10,'ok')
hold on; plot(FingersKinInfo.Finger(3).Es_MaxVel,10,'ob')

figure;
plot(Posq_F3(:,1))
hold on; plot(FingersKinInfo.Finger(3).Fs_MaxPos,-30,'og')
hold on; plot(FingersKinInfo.Finger(3).Es_MaxPos,-5,'or')
hold on; vline(FingersKinInfo.Finger(3).Fs_MaxVel,'g--')
hold on; vline(FingersKinInfo.Finger(3).Es_MaxVel,'r--')

FingersKinData.Finger(3).PurePos=Posq_F3(:,1);
FingersKinData.Finger(3).PureVel=PC_Vel(3).Movement;

%% Finger 4 position
half_window=37; %points; about 1 sec: for finding local max and local min in position data

% for finger 4 
Movement = 4;
Posq_F4 = PC_Pos(Movement).Movement;
[F4_n,F4_m]=size(Posq_F4);

% finding max and min time stamps
j=1;
jj=1;
for i=2*half_window+10:F4_n-(half_window+1)
    
    logic1=(Posq_F4(i-half_window:i+half_window,1)> Posq_F4(i,1));
    
    if sum(logic1)==2*half_window
        FingersKinInfo.Finger(4).Fs_MaxPos(j,1)=i;
        j=j+1;
    end
    
    logic2=(Posq_F4(i-half_window:i+half_window,1)< Posq_F4(i,1));
    
    if sum(logic2)==2*half_window
        FingersKinInfo.Finger(4).Es_MaxPos(jj,1)=i;
        jj=jj+1;
    end
    
    
end

% check on the plot
figure;
plot(Posq_F4(:,1))
hold on; plot(FingersKinInfo.Finger(4).Fs_MaxPos,-10,'og')
hold on; plot(FingersKinInfo.Finger(4).Es_MaxPos,-5,'or')

%% Finger 4 Velocity

Index_F4=sort([1; FingersKinInfo.Finger(4).Fs_MaxPos;...
    FingersKinInfo.Finger(4).Es_MaxPos; F4_n]);

AbsVel_F4 = abs(PC_Vel(4).Movement); 

j=1;
for i=1:2:(length(Index_F4)-1)
    Section=AbsVel_F4(Index_F4(i):(Index_F4(i+1)));
    [m_S n_S]=max(Section);
    FingersKinInfo.Finger(4).Fs_MaxVel(j,1)=(n_S+Index_F4(i)-1);  
    j=j+1;
end

j=1;
for i=2:2:(length(Index_F4)-1)
    Section=AbsVel_F4(Index_F4(i):(Index_F4(i+1)));
    [m_S n_S]=max(Section);
    FingersKinInfo.Finger(4).Es_MaxVel(j,1)=(n_S+Index_F4(i)-1);
    j=j+1;
end


% check on the plot
figure;
plot(AbsVel_F4)
hold on; plot(FingersKinInfo.Finger(4).Fs_MaxVel,10,'ok')
hold on; plot(FingersKinInfo.Finger(4).Es_MaxVel,10,'ob')

figure;
plot(Posq_F4(:,1))
hold on; plot(FingersKinInfo.Finger(4).Fs_MaxPos,-30,'og')
hold on; plot(FingersKinInfo.Finger(4).Es_MaxPos,-5,'or')
hold on; vline(FingersKinInfo.Finger(4).Fs_MaxVel,'g--')
hold on; vline(FingersKinInfo.Finger(4).Es_MaxVel,'r--')

FingersKinData.Finger(4).PurePos=Posq_F4(:,1);
FingersKinData.Finger(4).PureVel=PC_Vel(4).Movement;

%% Finger 5 position
half_window=37; %points; about 1 sec: for finding local max and local min in position data

% for finger 5 
Movement = 5;
Posq_F5 = PC_Pos(Movement).Movement;
[F5_n,F5_m]=size(Posq_F5);

% finding max and min time stamps
j=1;
jj=1;
for i=2*half_window+10:F5_n-(half_window+1)
    
    logic1=(Posq_F5(i-half_window:i+half_window,1)> Posq_F5(i,1));
    
    if sum(logic1)==2*half_window
        FingersKinInfo.Finger(5).Fs_MaxPos(j,1)=i;
        j=j+1;
    end
    
    logic2=(Posq_F5(i-half_window:i+half_window,1)< Posq_F5(i,1));
    
    if sum(logic2)==2*half_window
        FingersKinInfo.Finger(5).Es_MaxPos(jj,1)=i;
        jj=jj+1;
    end
    
    
end

% check on the plot
figure;
plot(Posq_F5(:,1))
hold on; plot(FingersKinInfo.Finger(5).Fs_MaxPos,-10,'og')
hold on; plot(FingersKinInfo.Finger(5).Es_MaxPos,-5,'or')

%% Finger 5 Velocity

Index_F5=sort([1; FingersKinInfo.Finger(5).Fs_MaxPos;...
    FingersKinInfo.Finger(5).Es_MaxPos; F5_n]);

AbsVel_F5 = abs(PC_Vel(5).Movement); 

j=1;
for i=1:2:(length(Index_F5)-1)
    Section=AbsVel_F5(Index_F5(i):(Index_F5(i+1)));
    [m_S n_S]=max(Section);
    FingersKinInfo.Finger(5).Fs_MaxVel(j,1)=(n_S+Index_F5(i)-1);  
    j=j+1;
end

j=1;
for i=2:2:(length(Index_F5)-1)
    Section=AbsVel_F5(Index_F5(i):(Index_F5(i+1)));
    [m_S n_S]=max(Section);
    FingersKinInfo.Finger(5).Es_MaxVel(j,1)=(n_S+Index_F5(i)-1);
    j=j+1;
end


% check on the plot
figure;
plot(AbsVel_F5)
hold on; plot(FingersKinInfo.Finger(5).Fs_MaxVel,10,'ok')
hold on; plot(FingersKinInfo.Finger(5).Es_MaxVel,10,'ob')

figure;
plot(Posq_F5(:,1))
hold on; plot(FingersKinInfo.Finger(5).Fs_MaxPos,-30,'og')
hold on; plot(FingersKinInfo.Finger(5).Es_MaxPos,-5,'or')
hold on; vline(FingersKinInfo.Finger(5).Fs_MaxVel,'g--')
hold on; vline(FingersKinInfo.Finger(5).Es_MaxVel,'r--')

FingersKinData.Finger(5).PurePos=-Posq_F5(:,1);
FingersKinData.Finger(5).PureVel=-PC_Vel(5).Movement;

%% Grasp 1 position
half_window=37; %points; about 1 sec: for finding local max and local min in position data

% for grasp 1 
Movement = 6;
Posq_F6 = PC_Pos(Movement).Movement;
[F6_n,F6_m]=size(Posq_F6);

% finding max and min time stamps
j=1;
jj=1;
for i=2*half_window+10:F6_n-(half_window+1)
    
    logic1=(Posq_F6(i-half_window:i+half_window,1)> Posq_F6(i,1));
    
    if sum(logic1)==2*half_window
        FingersKinInfo.Finger(6).Fs_MaxPos(j,1)=i;
        j=j+1;
    end
    
    logic2=(Posq_F6(i-half_window:i+half_window,1)< Posq_F6(i,1));
    
    if sum(logic2)==2*half_window
        FingersKinInfo.Finger(6).Es_MaxPos(jj,1)=i;
        jj=jj+1;
    end
    
    
end

% check on the plot
figure;
plot(Posq_F6(:,1))
hold on; plot(FingersKinInfo.Finger(6).Fs_MaxPos,-10,'og')
hold on; plot(FingersKinInfo.Finger(6).Es_MaxPos,-5,'or')

%% Grasp 1 Velocity

Index_F6=sort([1; FingersKinInfo.Finger(6).Fs_MaxPos;...
    FingersKinInfo.Finger(6).Es_MaxPos; F6_n]);

AbsVel_F6 = abs(PC_Vel(6).Movement); 

j=1;
for i=1:2:(length(Index_F6)-1)
    Section=AbsVel_F6(Index_F6(i):(Index_F6(i+1)));
    [m_S n_S]=max(Section);
    FingersKinInfo.Finger(6).Fs_MaxVel(j,1)=(n_S+Index_F6(i)-1);  
    j=j+1;
end

j=1;
for i=2:2:(length(Index_F6)-1)
    Section=AbsVel_F6(Index_F6(i):(Index_F6(i+1)));
    [m_S n_S]=max(Section);
    FingersKinInfo.Finger(6).Es_MaxVel(j,1)=(n_S+Index_F6(i)-1);
    j=j+1;
end


% check on the plot
figure;
plot(AbsVel_F6)
hold on; plot(FingersKinInfo.Finger(6).Fs_MaxVel,10,'ok')
hold on; plot(FingersKinInfo.Finger(6).Es_MaxVel,10,'ob')

figure;
plot(Posq_F6(:,1))
hold on; plot(FingersKinInfo.Finger(6).Fs_MaxPos,-30,'og')
hold on; plot(FingersKinInfo.Finger(6).Es_MaxPos,-5,'or')
hold on; vline(FingersKinInfo.Finger(6).Fs_MaxVel,'g--')
hold on; vline(FingersKinInfo.Finger(6).Es_MaxVel,'r--')

FingersKinData.Finger(6).PurePos=-Posq_F6(:,1);
FingersKinData.Finger(6).PureVel=-PC_Vel(6).Movement;


%% Grasp 2 position
half_window=37; %points; about 1 sec: for finding local max and local min in position data

% for grasp 2 
Movement = 7;
Posq_F7 = PC_Pos(Movement).Movement;
[F7_n,F7_m]=size(Posq_F7);

% finding max and min time stamps
j=1;
jj=1;
for i=2*half_window+10:F7_n-(half_window+1)
    
    logic1=(Posq_F7(i-half_window:i+half_window,1)> Posq_F7(i,1));
    
    if sum(logic1)==2*half_window
        FingersKinInfo.Finger(7).Fs_MaxPos(j,1)=i;
        j=j+1;
    end
    
    logic2=(Posq_F7(i-half_window:i+half_window,1)< Posq_F7(i,1));
    
    if sum(logic2)==2*half_window
        FingersKinInfo.Finger(7).Es_MaxPos(jj,1)=i;
        jj=jj+1;
    end
    
    
end

% check on the plot
figure;
plot(Posq_F7(:,1))
hold on; plot(FingersKinInfo.Finger(7).Fs_MaxPos,-10,'og')
hold on; plot(FingersKinInfo.Finger(7).Es_MaxPos,-5,'or')

%% Grasp 2 Velocity

Index_F7=sort([1; FingersKinInfo.Finger(7).Fs_MaxPos;...
    FingersKinInfo.Finger(7).Es_MaxPos; F7_n]);

AbsVel_F7 = abs(PC_Vel(7).Movement); 

j=1;
for i=1:2:(length(Index_F7)-1)
    Section=AbsVel_F7(Index_F7(i):(Index_F7(i+1)));
    [m_S n_S]=max(Section);
    FingersKinInfo.Finger(7).Fs_MaxVel(j,1)=(n_S+Index_F7(i)-1);  
    j=j+1;
end

j=1;
for i=2:2:(length(Index_F7)-1)
    Section=AbsVel_F7(Index_F7(i):(Index_F7(i+1)));
    [m_S n_S]=max(Section);
    FingersKinInfo.Finger(7).Es_MaxVel(j,1)=(n_S+Index_F7(i)-1);
    j=j+1;
end


% check on the plot
figure;
plot(AbsVel_F7)
hold on; plot(FingersKinInfo.Finger(7).Fs_MaxVel,10,'ok')
hold on; plot(FingersKinInfo.Finger(7).Es_MaxVel,10,'ob')

figure;
plot(Posq_F7(:,1))
hold on; plot(FingersKinInfo.Finger(7).Fs_MaxPos,-30,'og')
hold on; plot(FingersKinInfo.Finger(7).Es_MaxPos,-5,'or')
hold on; vline(FingersKinInfo.Finger(7).Fs_MaxVel,'g--')
hold on; vline(FingersKinInfo.Finger(7).Es_MaxVel,'r--')

FingersKinData.Finger(7).PurePos= Posq_F7(:,1);
FingersKinData.Finger(7).PureVel= PC_Vel(7).Movement;

%% Grasp 3 position
half_window=37; %points; about 1 sec: for finding local max and local min in position data

% for grasp 3 
Movement = 8;
Posq_F8 = PC_Pos(Movement).Movement;
[F8_n,F8_m]=size(Posq_F8);

% finding max and min time stamps
j=1;
jj=1;
for i=2*half_window:F8_n-(half_window+1)
    
    logic1=(Posq_F8(i-half_window:i+half_window,1)> Posq_F8(i,1));
    
    if sum(logic1)==2*half_window
        FingersKinInfo.Finger(8).Fs_MaxPos(j,1)=i;
        j=j+1;
    end
    
    logic2=(Posq_F8(i-half_window:i+half_window,1)< Posq_F8(i,1));
    
    if sum(logic2)==2*half_window
        FingersKinInfo.Finger(8).Es_MaxPos(jj,1)=i;
        jj=jj+1;
    end
    
    
end

% check on the plot
figure;
plot(Posq_F8(:,1))
hold on; plot(FingersKinInfo.Finger(8).Fs_MaxPos,-10,'og')
hold on; plot(FingersKinInfo.Finger(8).Es_MaxPos,-5,'or')

%% Grasp 3 Velocity

Index_F8=sort([1; FingersKinInfo.Finger(8).Fs_MaxPos;...
    FingersKinInfo.Finger(8).Es_MaxPos; F8_n]);

AbsVel_F8 = abs(PC_Vel(8).Movement); 

j=1;
for i=1:2:(length(Index_F8)-1)
    Section=AbsVel_F8(Index_F8(i):(Index_F8(i+1)));
    [m_S n_S]=max(Section);
    FingersKinInfo.Finger(8).Fs_MaxVel(j,1)=(n_S+Index_F8(i)-1);  
    j=j+1;
end

j=1;
for i=2:2:(length(Index_F8)-1)
    Section=AbsVel_F8(Index_F8(i):(Index_F8(i+1)));
    [m_S n_S]=max(Section);
    FingersKinInfo.Finger(8).Es_MaxVel(j,1)=(n_S+Index_F8(i)-1);
    j=j+1;
end

% check on the plot
figure;
plot(AbsVel_F8)
hold on; plot(FingersKinInfo.Finger(8).Fs_MaxVel,10,'ok')
hold on; plot(FingersKinInfo.Finger(8).Es_MaxVel,10,'ob')

figure;
plot(Posq_F8(:,1))
hold on; plot(FingersKinInfo.Finger(8).Fs_MaxPos,-30,'og')
hold on; plot(FingersKinInfo.Finger(8).Es_MaxPos,-5,'or')
hold on; vline(FingersKinInfo.Finger(8).Fs_MaxVel,'g--')
hold on; vline(FingersKinInfo.Finger(8).Es_MaxVel,'r--')

FingersKinData.Finger(8).PurePos= Posq_F8(:,1);
FingersKinData.Finger(8).PureVel= PC_Vel(8).Movement;


%% process of neural data

% calculating HG-LFO ;automatically will be hilbert value
% only the first 256 channels

for Fi=1:8
    Signal=hg{1,Fi};
    
    % for hg-lfo
    F_Envelope = Signal(:,65:320); 
    [b,a]=butter(3,[.5, 4]/(Fs/2));
    HG_mean = mean(F_Envelope,1);
    hglfo_hilbert{1,Fi}=filtfilt(b,a,F_Envelope)+repmat(HG_mean,length(F_Envelope),1);
    
end

% calculating LFO ;automatically will be hilbert value
for Fi=1:8
    Signal=delta{1,Fi};
    lfo_hilbert{1,Fi} = Signal(:,65:320);
 
end

% zscore across all fingers and grasp 
% for all finger
Cat_LFO = [];
Cat_HG = [];
for Fi = 1:8
    Cat_LFO = [Cat_LFO;  lfo_hilbert{1,Fi}];
    Cat_HG = [Cat_HG; hglfo_hilbert{1,Fi}];
end
LFO_Hilbert_all = zscore(Cat_LFO);
HG_Hilbert_all = zscore(Cat_HG);

% return data to the same structures
Finger_Sizes = [1, length(delta{1,1}),length(delta{1,2}),length(delta{1,3}),length(delta{1,4}),length(delta{1,5})...
    length(delta{1,6}), length(delta{1,7}), length(delta{1,8})];

for Fi = 1:8
    
    lfo_hilbert_zscore{1,Fi} = LFO_Hilbert_all(sum(Finger_Sizes(1:Fi)):sum(Finger_Sizes(1:Fi+1))-1,:);
    hglfo_hilbert_zscore{1,Fi} = HG_Hilbert_all(sum(Finger_Sizes(1:Fi)):sum(Finger_Sizes(1:Fi+1))-1,:);
    
end 


%% find the neural data (window of trials) for max vel for flexion and extension

% finding the trial periods for each finger and each grasp, save for dPCA 
Win = 25;

% for lfo for flexion

for Fi = 1:8
        
    for i = 1:length(FingersKinInfo.Finger(Fi).Fs_MaxVel)-1
       
        index = FingersKinInfo.Finger(Fi).Fs_MaxVel(i,1);
        
        Trial_Data = lfo_hilbert_zscore{1, Fi}(index-Win:index+Win, 1:256);
       
        ECoG_lfo(Fi).finger(i).Trials = Trial_Data;
         
    end
      
end 

% for hglfo for flexion

for Fi = 1:8
        
    for i = 1:length(FingersKinInfo.Finger(Fi).Fs_MaxVel)-1
       
        index = FingersKinInfo.Finger(Fi).Fs_MaxVel(i,1);
        
        Trial_Data = hglfo_hilbert_zscore{1, Fi}(index-Win:index+Win, 1:256);
       
         ECoG_hglfo(Fi).finger(i).Trials = Trial_Data;
         
    end
      
end 

%% save related data

save('/media/reza/WindowsDrive/ECoGLeapMotion/ResultsGroupAnalysis/github_Branch_V3/Subject2_NikMethod_25Hz.mat',...
    'ECoG_lfo','ECoG_hglfo');
