clear;

addpath(genpath('../3rdParty'));
addpath(genpath('../matlab'));
dataPath = '~/research/data/Hopkins155_AdditionalSequences_MissingData'; 

file = dir(dataPath);
ii = 0;
% rng(0);
% for i = 1:length(file)
for i = 9
	if( (file(i).isdir == 1) && ~strcmp(file(i).name,'.') && ~strcmp(file(i).name,'..') )
		filePath = file(i).name;
% 		eval(['cd ' filepath]);
        f = dir(fullfile(dataPath, filePath));
		
        foundValidData = false;
        for j = 1:length(f)
            if( ~isempty(strfind(f(j).name,'_truth.mat')) )
                ind = j;
                foundValidData = true;
                load(fullfile(dataPath, filePath, f(ind).name));
                break
            end
        end
		
		if(foundValidData)
            
            [s, ind] = sort(s); x = x(:, ind, :); y = y(:, ind, :); m = m(ind, :);
            theta = pi / 4; ratio = 0.5;
            [x, y, camID] = dataRotation(x, y, K, theta, ratio);
            [camID, idx] = sort(camID); x = x(:, idx, :); y = y(:, idx, :); s = s(idx); m = m(idx, :);
            
			N = size(x,2);
			F = size(x,3);
			D = 2*F;	
						
            X = reshape(permute(x(1:2,:,:),[1 3 2]),D,N);
            

            
            % remove short trajectories
            minLen = 5;
            indShort = getShortTrajIndex(m, minLen);
            X(:, indShort) = [];
            m(indShort, :) = [];
            s(indShort) = [];
            camID(indShort) = [];
            indX = kron(m,[1,1])';		
			
			X = X.*indX;			
			
            rng('default');
%             [missrate, grp, bestRank, minNcutValue] = RSIM_Incomplete(X',indX', s, 4, 1);
% 			[missrate, grp, bestRank, minNcutValue] = RSIMJBLD_Incomplete(X',indX',s,4,1);
            [missrate, grp, bestRank,W, index] = imprvRSIM_Incomplete(X', indX', s, 4, 1, camID);
%             [missrate, grp, bestRank,W, index] = imprvRSIM_JBLD_Incomplete(X', indX', s, 4, 1, camID);
           

			ii = ii+1;
			Missrate(ii) = missrate;
			disp([filePath ': ' num2str(100*Missrate(ii)) '%' ', nMotions: ' num2str(max(s)) ', seq: ' num2str(ii)]);
%             clear x y s m
		end
	end
end

avgtol = mean(Missrate);
medtol = median(Missrate);
maxtol = max(Missrate);
stdtol = std(Missrate);

disp('Results on Hopkins Missing data')
disp(['Mean: ' num2str(100*avgtol) '%' ', median: ' num2str(100*medtol) '%;'...
      ', max: ' num2str(100*maxtol) '%;' ', std: ' num2str(100*stdtol) '%;']);
