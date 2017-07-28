


% Motion segmentation with Robust Shape Interaction Matrix Method with JBLD

clear; close all
addpath(genpath('../3rdParty'));
addpath(genpath('../matlab'));

% load('toyTrajectories_tranANDrot')
rotFreq1 = 360;
rotFreq2 = 75;
vel1freq = 9/2;
vel2freq = 11/2;
[imPointsc1, imPointsc2sync, imPointsc2, s1,s2,s3] = ...
    toyTrajectories(rotFreq1, rotFreq2, vel1freq, vel2freq);


% file = dir(dataPath);
ii = 0;
ii2 = 0;
ii3 = 0;
tic

for i = 1:1%length(file)
    
    y1= imPointsc1;
    x1 = [(y1(1,:,:) / (max(max(y1(1,:,:)))/2))-1;(y1(2,:,:) / (max(max(y1(2,:,:)))/2))-1;y1(3,:,:)];
    y2= imPointsc2;
    x2 = [(y2(1,:,:) / (max(max(y2(1,:,:)))/2))-1;(y2(2,:,:) / (max(max(y2(2,:,:)))/2))-1;y2(3,:,:)];
    x= cat(2,y1(:,1:4,:), y2(:,1:4,:), y1(:,5:8,:) ,y2(:,5:8,:));
    s = [s1(1:4) ; s2(1:4);s1(5:8) ; s2(5:8)];
    
%     x=x1;
%     s = s2;
%     
    N = size(x,2);
    F = size(x,3);
    D = 3*F;
    
    X = reshape(permute(x(1:3,:,:),[1 3 2]),D,N);	% note here the all-one rows are also included
    
    % 			[missrate, grp, bestRank, minNcutValue,W] = RSIM(X, s);
    [missrate, grp, bestRank, minNcutValue,W] = RSIM_JBLD(X, s);
    
    ii = ii+1;
    Missrate(ii) = missrate;
    disp([num2str(100*Missrate(ii)) '%, dim:' num2str(bestRank) ', nMotions: ' num2str(max(s)) ', seq: ' num2str(ii)]);
    if(max(s)==2)
        ii2 = ii2+1;
        Missrate2(ii2) = Missrate(ii);
    else
        ii3 = ii3+1;
        Missrate3(ii3) = Missrate(ii);
    end
    
end

time = toc;
avgtime = time/ii

avgtol = mean(Missrate);
medtol = median(Missrate);
avgtwo = mean(Missrate2);
medtwo = median(Missrate2);
% avgthree = mean(Missrate3);
% medthree = median(Missrate3);

disp('Results on Hopkins155')
disp(['Mean of all: ' num2str(100*avgtol) '%' ', median of all: ' num2str(100*medtol) '%;']);
disp(['Mean of two: ' num2str(100*avgtwo) '%' ', median of two: ' num2str(100*medtwo) '%;']);
% disp(['Mean of three: ' num2str(100*avgthree) '%' ', median of three: ' num2str(100*medthree) '%.']);