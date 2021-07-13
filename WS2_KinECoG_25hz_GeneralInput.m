% analyzing the Kinematics data and ECoG data. 
% in this version; save:
% 1- frequency of recording  
% 2- the index for pos max and vel max for later sepearting flexion and extension
% 3- the brain areas; layout for all good channels 
% 4- related processed ECoG signals

close all;
clear all;
clc;

%% using Nik preprocessed data: Generating ECoG data per trial
load('/media/reza/WindowsDrive/ECoGLeapMotion/DataPatientTwo/AllData/ecog_kinematics_25Hz.mat')

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


% finding the ECoG trial periods for each finger 

% for lfo
for Fi = 1:5
    
    time_stamp = tstart{1, Fi};
    
    for i = 1:length(time_stamp)-2
        Trial_Interp = zeros(75,256);
        
        index = find(time{1, Fi}>time_stamp(i) & time{1, Fi}<time_stamp(i+1));
        
        Trial_Data = lfo_hilbert{1, Fi}(index, 1:256);
       
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
       
       ECoGKin_lfo(Fi).FingerECoG(i).Trials = Trial_Interp;
         
    end 
    
end 

% for hglfo
for Fi = 1:5
    
    time_stamp = tstart{1, Fi};
    
    for i = 1:length(time_stamp)-2
        Trial_Interp = zeros(75,256);
        
        index = find(time{1, Fi}>time_stamp(i) & time{1, Fi}<time_stamp(i+1));
        
        Trial_Data = hglfo_hilbert{1, Fi}(index, 1:256);
       
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
       
       ECoGKin_hglfo(Fi).FingerECoG(i).Trials = Trial_Interp;
         
    end 
    
end 


%% using Nik preprocessed data: Generating finger Kinematics data per trial

for Fi = 1:5
    
    % using PCA for kinematics data
    Data_pos = [];
    Data_vel = [];
    for Finger = 1:5
        
        Data_pos = [Data_pos fingers{1, Fi}(Finger).spos];
        Data_vel = [Data_vel fingers{1, Fi}(Finger).svel];
        
    end
    
    [Pos_coeff,Pos_score,Pos_latent,Pos_tsquared,Pos_explained] = pca(Data_pos);
    [Vel_coeff,Vel_score,Vel_latent,Vel_tsquared,Vel_explained] = pca(Data_vel);
    
    time_stamp = tstart{1, Fi};
    
    % for Pos_score
    for i = 1:length(time_stamp)-2
        Trial_Interp = zeros(75,15);
        
        index = find(time{1, Fi}>time_stamp(i) & time{1, Fi}<time_stamp(i+1));
        
        Trial_Data = Pos_score(index, :);
        
        if length(index) == 75
            
            Trial_Interp = Trial_Data;
            
        else
            xq = 0:(length(index)/75):length(index)-1;
            if length(xq) == 75
                xq = xq;
            else
                xq = [xq,length(index)-1];
            end
            
            for j=1:15
                
                Trial_Interp(:,j) = interp1(0:1:length(index)-1,Trial_Data(:,j) ,xq);
            end
            
        end
        
        ECoGKin_hglfo(Fi).FingerPCAPos(i).Trials = Trial_Interp;
        
    end
    
    % for Vel_score
    for i = 1:length(time_stamp)-2
        Trial_Interp = zeros(75,15);
        
        index = find(time{1, Fi}>time_stamp(i) & time{1, Fi}<time_stamp(i+1));
        
        Trial_Data = Vel_score(index, :);
        
        if length(index) == 75
            
            Trial_Interp = Trial_Data;
            
        else
            xq = 0:(length(index)/75):length(index)-1;
            if length(xq) == 75
                xq = xq;
            else
                xq = [xq,length(index)-1];
            end
            
            for j=1:15
                
                Trial_Interp(:,j) = interp1(0:1:length(index)-1,Trial_Data(:,j) ,xq);
            end
            
        end
        
        ECoGKin_hglfo(Fi).FingerPCAVel(i).Trials = Trial_Interp;
        
    end
    
    
end



%% save related data
% Fs, FingersKinData, FingersKinInfo, Selected_Chs,    

save('/media/reza/WindowsDrive/ECoGLeapMotion/ResultsGroupAnalysis/github_Branch_V3/WS2_KinECoG_25hz.mat',...
    'ECoGKin_lfo','ECoGKin_hglfo');


