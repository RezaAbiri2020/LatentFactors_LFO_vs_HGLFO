% analyzing the Kinematics data and ECoG data. 
% in this version; save:
% 1- frequency of recording  
% 2- the index for pos max and vel max for later sepearting flexion and extension
% 3- the brain areas; layout for all good channels 
% 4- related processed ECoG signals

close all;
clear all;
clc;

%% using Nik preprocessed data

load('/media/reza/WindowsDrive/ECoGLeapMotion/DataPatientOne/AllData/ecog_kinematics_25Hz.mat')

% calculating HG-LFO ;automatically will be hilbert value
Fs = 25;

for Fi=1:5
    Signal=hg{1,Fi};
    
    % for hg-lfo
    %[F_Envelope,Lower]=envelope(Signal);
    F_Envelope = Signal; 
    [b,a]=butter(3,[.5, 4]/(Fs/2));
    HG_mean = mean(F_Envelope,1);
    hglfo_hilbert{1,Fi}=filtfilt(b,a,F_Envelope)+repmat(HG_mean,length(F_Envelope),1);
    
end

% hilbert of delta
for Fi=1:5
    Signal=delta{1,Fi};
    lfo_hilbert{1,Fi} = abs(hilbert(Signal));
 
end

% zscore across all fingers
% for all finger
Cat_LFO = [];
Cat_HG = [];
for Fi = 1:5
    Cat_LFO = [Cat_LFO;  lfo_hilbert{1,Fi}];
    Cat_HG = [Cat_HG; hglfo_hilbert{1,Fi}];
end
LFO_Hilbert_all = zscore(Cat_LFO);
HG_Hilbert_all = zscore(Cat_HG);

% return data to the same structures
Finger_Sizes = [1, length(delta{1,1}),length(delta{1,2}),length(delta{1,3}),length(delta{1,4}),length(delta{1,5})];

for Fi = 1:5
    
    lfo_hilbert_zscore{1,Fi} = LFO_Hilbert_all(sum(Finger_Sizes(1:Fi)):sum(Finger_Sizes(1:Fi+1))-1,:);
    hglfo_hilbert_zscore{1,Fi} = HG_Hilbert_all(sum(Finger_Sizes(1:Fi)):sum(Finger_Sizes(1:Fi+1))-1,:);
    
end 




% finding the trial periods for each finger, save for dPCA 

% for lfo
for Fi = 1:5
    
    time_stamp = tstart{1, Fi};
    
    for i = 1:length(time_stamp)-2
        Trial_Interp = zeros(75,256);
        
        index = find(time{1, Fi}>time_stamp(i) & time{1, Fi}<time_stamp(i+1));
        
        Trial_Data = lfo_hilbert_zscore{1, Fi}(index, 1:256);
       
       if length(index) == 75
           
           Trial_Interp = Trial_Data;
           
       else
           xq = 0:(length(index)/75):length(index)-1;
           if length(xq) == 75
               xq = xq;
           else 
               xq = [xq,length(index)-1];
           end 
           
           for j=1:256  
               
               Trial_Interp(:,j) = interp1(0:1:length(index)-1,Trial_Data(:,j) ,xq);          
           end 
           
       end
       
       ECoG_lfo(Fi).finger(i).Trials = Trial_Interp;
         
    end 
    
end 

% for hglfo
for Fi = 1:5
    
    time_stamp = tstart{1, Fi};
    
    for i = 1:length(time_stamp)-2
        Trial_Interp = zeros(75,256);
        
        index = find(time{1, Fi}>time_stamp(i) & time{1, Fi}<time_stamp(i+1));
        
        Trial_Data = hglfo_hilbert_zscore{1, Fi}(index, 1:256);
       
       if length(index) == 75
           
           Trial_Interp = Trial_Data;
           
       else
           xq = 0:(length(index)/75):length(index)-1;
           if length(xq) == 75
               xq = xq;
           else 
               xq = [xq,length(index)-1];
           end 
           
           for j=1:256  
               
               Trial_Interp(:,j) = interp1(0:1:length(index)-1,Trial_Data(:,j) ,xq);          
           end 
           
       end
       
       ECoG_hglfo(Fi).finger(i).Trials = Trial_Interp;
         
    end 
    
end 



%% save related data
% Fs, FingersKinData, FingersKinInfo, Selected_Chs,    

%save('/media/reza/WindowsDrive/ECoGLeapMotion/ResultsGroupAnalysis/github_Branch_V3/Subject1_NikMethod.mat',...
%    'ECoG_lfo','ECoG_hglfo');

