function [trajectory] = ball2trajectory(D)
numFrames = size(D,1);
trajectory = nan(numFrames+1,3);
trajectory(1,:) = [pi/2,0,0];
mmConst  = 8.32;
radConst = mmConst/3;
D(:,1) = D(:,1) *  mmConst; % x values
D(:,2) = D(:,2) *  mmConst; % y values
D(:,3) = D(:,3) * radConst; % z values
for iFrame = 1:numFrames
    currentPos = trajectory(iFrame,:);
    currentPos(1) = currentPos(1) + (D(iFrame,3) / 50) * -1;
    currentPos(2) = currentPos(2) +...
                    (D(iFrame,1)/50) * cos(currentPos(1)) -...
                    (D(iFrame,2)/50) * sin(currentPos(1));
    
    currentPos(3) = currentPos(3) +...
                    (D(iFrame,1)/50) * sin(currentPos(1)) -...
                    (D(iFrame,2)/50) * cos(currentPos(1));
    trajectory(iFrame+1,:) = currentPos;
end