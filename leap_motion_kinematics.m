function [time,hand,fingers] = leap_motion_kinematics(data)

% sort fingers
for frame = 1:length(data),
    if isempty(data(frame).pointables),
        continue;
    else,
        fingers = data(frame).pointables; 
        finger_id = [fingers(1).id fingers(2).id fingers(3).id fingers(4).id fingers(5).id];
        if issorted(finger_id), % fingers are in order
            continue;
        else, % fix the order of the fingers
            [~, sorted_index] = sort(finger_id);
            
            finger1 = data(frame).pointables(sorted_index(1));            
            finger2 = data(frame).pointables(sorted_index(2));
            finger3 = data(frame).pointables(sorted_index(3));
            finger4 = data(frame).pointables(sorted_index(4));
            finger5 = data(frame).pointables(sorted_index(5));
            
            data(frame).pointables(1) = finger1;
            data(frame).pointables(2) = finger2;
            data(frame).pointables(3) = finger3;
            data(frame).pointables(4) = finger4; 
            data(frame).pointables(5) = finger5;
            
        end
    end  
end
clear fingers

% remove frames without 1 hand and 5 fingers
for i=1:length(data),
    hand_idx(i) = length(data(i).hands)==1;
    if hand_idx(i),
        fingers_idx(i) = length(data(i).pointables)==5;
    else
        fingers_idx(i) = false;
    end
end
idx = hand_idx & fingers_idx;
data = data(idx);

% remove repeated frames
idx = [];
prev_id = 0;
for i=1:length(data),
    id = data(i).id;
    if id==prev_id,
        idx(i) = false;
    else
        idx(i) = true;
    end
    prev_id = id;
end
data = data(idx==1);

% hand kinematics
tmp = cat(1,data.hands);
hand.palm_pos = cat(1,tmp.palm_position);
hand.palm_vel = cat(1,tmp.palm_velocity);
hand.palm_normal = cat(1,tmp.palm_normal);
hand.basis = cat(3,tmp.basis);

% finger kinematics
tmp = cat(1,data.pointables);
% fingers(1).pos = cat(1,tmp(:,1).position);
% fingers(1).dir = cat(1,tmp(:,1).direction);
% fingers(2).pos = cat(1,tmp(:,2).position);
% fingers(2).dir = cat(1,tmp(:,2).direction);
% fingers(3).pos = cat(1,tmp(:,3).position);
% fingers(3).dir = cat(1,tmp(:,3).direction);
% fingers(4).pos = cat(1,tmp(:,4).position);
% fingers(4).dir = cat(1,tmp(:,4).direction);
% fingers(5).pos = cat(1,tmp(:,5).position);
% fingers(5).dir = cat(1,tmp(:,5).direction);

% finger positions in hand basis
for frame=1:length(data),
    for f=1:5,
        fingers(f).pos(frame,:) = (tmp(frame,f).position - hand.palm_pos(frame,:)) * hand.basis(:,:,frame);
        fingers(f).vel(frame,:) = tmp(frame,f).velocity * hand.basis(:,:,frame);
        fingers(f).dir(frame,:) = tmp(frame,f).direction * hand.basis(:,:,frame);
    end
end

% timestamps
time = cat(1,data.timestamp) * 1e-6; % convert from usecs to secs
time = time - time(1); % relative time

% remove avg from each pos and vel

