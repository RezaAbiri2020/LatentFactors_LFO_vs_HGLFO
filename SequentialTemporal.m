
% for this commit:
% looking to the ERP based on max pos or vel pos
% looking to the first half flexion (or extension) trials and compared to the second half

clear all
close all 
clc

%% load the data and indexing

% signals:
load('E:\ECoGLeapMotion\DataPatientTwo\github_Branch_V1\LFO_signals.mat')
load('E:\ECoGLeapMotion\DataPatientTwo\github_Branch_V1\HG_Direct_LFO_Signals.mat')
load('E:\ECoGLeapMotion\DataPatientTwo\github_Branch_V1\HG_Avg_LFO_Signals.mat')

% using indexing for flexsion/extension
load('E:\ECoGLeapMotion\DataPatientTwo\github_Branch_V1\FingersKinIndexes.mat')


%% ERP of max pos calculation for flexion

%choose a window
win=256; %sample points

% hilbert; organzie and do averaging for data 
for Fi=1:5
    Indexes=FingersKinIndexes.Finger(Fi).Fs_MaxPos;
    MaxPos_Fs_LFO=zeros(win,256,length(Indexes));
    MaxPos_Fs_DirectHG=zeros(win,256,length(Indexes));
    MaxPos_Fs_AvgHG=zeros(win,256,length(Indexes));
    
    Amp_LFO=abs(hilbert(LFO_signals(Fi).Filtered));
    Amp_Direct=abs(hilbert(HG_Direct_LFO_Signals(Fi).DeltaofEnv));
    Amp_Avg=abs(hilbert(HG_Avg_LFO_Signals(Fi).DeltaofEnv));
    
    for i=1:length(Indexes)
        MaxPos_Fs_LFO(:,:,i)=Amp_LFO((Indexes(i)-win+1):Indexes(i),:);
        MaxPos_Fs_DirectHG(:,:,i)=Amp_Direct((Indexes(i)-win+1):Indexes(i),:);
        MaxPos_Fs_AvgHG(:,:,i)=Amp_Avg((Indexes(i)-win+1):Indexes(i),:);
    end
    
    MaxPos_Fs_LFO_ERP1(Fi,:,:)=squeeze(mean(MaxPos_Fs_LFO(:,:,1:10),3));
    MaxPos_Fs_DirectHG_ERP1(Fi,:,:)=squeeze(mean(MaxPos_Fs_DirectHG(:,:,1:10),3));
    MaxPos_Fs_AvgHG_ERP1(Fi,:,:)=squeeze(mean(MaxPos_Fs_AvgHG(:,:,1:10),3));
    
    MaxPos_Fs_LFO_ERP2(Fi,:,:)=squeeze(mean(MaxPos_Fs_LFO(:,:,11:end),3));
    MaxPos_Fs_DirectHG_ERP2(Fi,:,:)=squeeze(mean(MaxPos_Fs_DirectHG(:,:,11:end),3));
    MaxPos_Fs_AvgHG_ERP2(Fi,:,:)=squeeze(mean(MaxPos_Fs_AvgHG(:,:,11:end),3));
        
end

%% ERP of max pos calculation for extension

%choose a window
win=256; %sample points

% hilbert; organzie and do averaging for data 
for Fi=1:5
    Indexes=FingersKinIndexes.Finger(Fi).Es_MaxPos;
    MaxPos_Es_LFO=zeros(win,256,length(Indexes));
    MaxPos_Es_DirectHG=zeros(win,256,length(Indexes));
    MaxPos_Es_AvgHG=zeros(win,256,length(Indexes));
    
    Amp_LFO=abs(hilbert(LFO_signals(Fi).Filtered));
    Amp_Direct=abs(hilbert(HG_Direct_LFO_Signals(Fi).DeltaofEnv));
    Amp_Avg=abs(hilbert(HG_Avg_LFO_Signals(Fi).DeltaofEnv));
    
    for i=1:length(Indexes)
        MaxPos_Es_LFO(:,:,i)=Amp_LFO((Indexes(i)-win+1):Indexes(i),:);
        MaxPos_Es_DirectHG(:,:,i)=Amp_Direct((Indexes(i)-win+1):Indexes(i),:);
        MaxPos_Es_AvgHG(:,:,i)=Amp_Avg((Indexes(i)-win+1):Indexes(i),:);
    end
    
    MaxPos_Es_LFO_ERP1(Fi,:,:)=squeeze(mean(MaxPos_Es_LFO(:,:,1:10),3));
    MaxPos_Es_DirectHG_ERP1(Fi,:,:)=squeeze(mean(MaxPos_Es_DirectHG(:,:,1:10),3));
    MaxPos_Es_AvgHG_ERP1(Fi,:,:)=squeeze(mean(MaxPos_Es_AvgHG(:,:,1:10),3));
    
    MaxPos_Es_LFO_ERP2(Fi,:,:)=squeeze(mean(MaxPos_Es_LFO(:,:,11:end),3));
    MaxPos_Es_DirectHG_ERP2(Fi,:,:)=squeeze(mean(MaxPos_Es_DirectHG(:,:,11:end),3));
    MaxPos_Es_AvgHG_ERP2(Fi,:,:)=squeeze(mean(MaxPos_Es_AvgHG(:,:,11:end),3));
        
end

%% ERP of max vel calculation for flexion

%choose a window
win=250; %sample points

% hilbert; organzie and do averaging for data 
for Fi=1:5
    Indexes=FingersKinIndexes.Finger(Fi).Fs_MaxVel;
    MaxVel_Fs_LFO=zeros(2*win,256,length(Indexes));
    MaxVel_Fs_DirectHG=zeros(2*win,256,length(Indexes));
    MaxVel_Fs_AvgHG=zeros(2*win,256,length(Indexes));
    
    Amp_LFO=abs(hilbert(LFO_signals(Fi).Filtered));
    Amp_Direct=abs(hilbert(HG_Direct_LFO_Signals(Fi).DeltaofEnv));
    Amp_Avg=abs(hilbert(HG_Avg_LFO_Signals(Fi).DeltaofEnv));
    
    for i=1:length(Indexes)
        MaxVel_Fs_LFO(:,:,i)=Amp_LFO((Indexes(i)-win+1):(Indexes(i)+win),:);
        MaxVel_Fs_DirectHG(:,:,i)=Amp_Direct((Indexes(i)-win+1):(Indexes(i)+win),:);
        MaxVel_Fs_AvgHG(:,:,i)=Amp_Avg((Indexes(i)-win+1):(Indexes(i)+win),:);
    end
    
    MaxVel_Fs_LFO_ERP1(Fi,:,:)=squeeze(mean(MaxVel_Fs_LFO(:,:,1:10),3));
    MaxVel_Fs_DirectHG_ERP1(Fi,:,:)=squeeze(mean(MaxVel_Fs_DirectHG(:,:,1:10),3));
    MaxVel_Fs_AvgHG_ERP1(Fi,:,:)=squeeze(mean(MaxVel_Fs_AvgHG(:,:,1:10),3));
    
    MaxVel_Fs_LFO_ERP2(Fi,:,:)=squeeze(mean(MaxVel_Fs_LFO(:,:,11:end),3));
    MaxVel_Fs_DirectHG_ERP2(Fi,:,:)=squeeze(mean(MaxVel_Fs_DirectHG(:,:,11:end),3));
    MaxVel_Fs_AvgHG_ERP2(Fi,:,:)=squeeze(mean(MaxVel_Fs_AvgHG(:,:,11:end),3));
        
end

%% ERP of max vel calculation for extension

%choose a window
win=250; %sample points

% hilbert; organzie and do averaging for data 
for Fi=1:5
    Indexes=FingersKinIndexes.Finger(Fi).Es_MaxVel;
    MaxVel_Es_LFO=zeros(2*win,256,length(Indexes));
    MaxVel_Es_DirectHG=zeros(2*win,256,length(Indexes));
    MaxVel_Es_AvgHG=zeros(2*win,256,length(Indexes));
    
    Amp_LFO=abs(hilbert(LFO_signals(Fi).Filtered));
    Amp_Direct=abs(hilbert(HG_Direct_LFO_Signals(Fi).DeltaofEnv));
    Amp_Avg=abs(hilbert(HG_Avg_LFO_Signals(Fi).DeltaofEnv));
    
    for i=1:length(Indexes)
        MaxVel_Es_LFO(:,:,i)=Amp_LFO((Indexes(i)-win+1):(Indexes(i)+win),:);
        MaxVel_Es_DirectHG(:,:,i)=Amp_Direct((Indexes(i)-win+1):(Indexes(i)+win),:);
        MaxVel_Es_AvgHG(:,:,i)=Amp_Avg((Indexes(i)-win+1):(Indexes(i)+win),:);
    end
    
    MaxVel_Es_LFO_ERP1(Fi,:,:)=squeeze(mean(MaxVel_Es_LFO(:,:,1:10),3));
    MaxVel_Es_DirectHG_ERP1(Fi,:,:)=squeeze(mean(MaxVel_Es_DirectHG(:,:,1:10),3));
    MaxVel_Es_AvgHG_ERP1(Fi,:,:)=squeeze(mean(MaxVel_Es_AvgHG(:,:,1:10),3));
    
    MaxVel_Es_LFO_ERP2(Fi,:,:)=squeeze(mean(MaxVel_Es_LFO(:,:,11:end),3));
    MaxVel_Es_DirectHG_ERP2(Fi,:,:)=squeeze(mean(MaxVel_Es_DirectHG(:,:,11:end),3));
    MaxVel_Es_AvgHG_ERP2(Fi,:,:)=squeeze(mean(MaxVel_Es_AvgHG(:,:,11:end),3));
        
end

%% fiding high peak for ERP1 and ERP2 and the ordering of peaks and plot the signals

% choose your data; it can be from ERP1 or ERP2
% choose a finger
Fi=1;
%ERP=squeeze(MaxPos_Fs_LFO_ERP2(Fi,:,:));
%ERP=squeeze(MaxPos_Fs_DirectHG_ERP2(Fi,:,:));

ERP=squeeze(MaxVel_Fs_LFO_ERP1(Fi,:,:));
%ERP=squeeze(MaxVel_Fs_DirectHG_ERP1(Fi,:,:));


data1=ERP;
data2=zeros(size(data1));
for i=1:256
    data2(:,i)=normalize(data1(:,i),'range');
end

data3=data2';

% figure;
% set(gcf, 'Position', [300, 100, 700, 850]);
% suptitle(['All Ch; MaxPos Fs LFO ERP2; Finger: ',num2str(Fi)]);
% imagesc(data3)
% colorbar

Rows=[];
for i=1:2*win
    [m,n]=max(data3(:,i));
    if isempty(find(Rows==n))
        Rows=[Rows;n];
    end
    
end

OrderedData=data3(Rows,:);
% figure;
% set(gcf, 'Position', [100, 100, 700, 850]);
% suptitle(['Ordered Ch; MaxPos Fs LFO ERP1; Finger: ',num2str(Fi)]);
% for i=1:10:100
%     subplot(10,1,ceil(i/10))
%     plot(OrderedData(i,:));
%     title(['Ch: ',num2str(Rows(i))])
%     
% end

figure;
set(gcf, 'Position', [300, 100, 700,850]);
suptitle(['Ordered Ch; MaxVel Fs LFO ERP2; Finger: ',num2str(Fi)]);
imagesc(OrderedData)
colorbar
yticks('')
yticks(1:3:length(Rows))
ytickangle(0)
yticklabels(num2str(Rows(1:3:end)))

% for real signals

% data4=ERP';
% OrderedSignal=data4(Rows,:);
% figure;
% set(gcf, 'Position', [300, 100, 700, 850]);
% suptitle(['Ordered Signal; MaxPos Fs LFO ERP2; Finger: ',num2str(Fi)]);
% for i=1:20:400
%     subplot(10,1,ceil(i/20))
%     plot(OrderedSignal(i,:));
%     title(['Ch: ',num2str(Rows(i))])
%     
% end
% 
% figure;
% set(gcf, 'Position', [300, 100, 700, 850]);
% suptitle(['Ordered Ch; MaxPos Fs LFO ERP2; Finger: ',num2str(Fi)]);
% imagesc(OrderedSignal)
% colorbar

%% map delta channel to HG

data1_new=HG_ERPHilbert(Fi).FingerF2E;
data2_new=zeros(size(data1));
    for i=1:256
        data2_new(:,i)=normalize(data1_new(:,i),'range');
    end
    
data2_new=data2_new';

    
OrderedData1_new=data2_new(Rows,:);
figure;
set(gcf, 'Position', [300, 100, 700, 850]);
%suptitle(['Ordered Ch; HG-LFO for Finger: ',num2str(Fi)]);
suptitle(['Ordered Ch of HG-LFO; projection of LFO Ch for Finger: ',num2str(Fi)]);
imagesc(OrderedData1_new)
colorbar
yticks('')
yticks(1:3:length(Rows))
ytickangle(0)
yticklabels(num2str(Rows(1:3:end)))
