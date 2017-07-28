
addpath(genpath('../3rdParty'));
addpath(genpath('../matlab'));
dataPath = '~/research/data/Hopkins155_AdditionalSequences_MissingData'; 

file = dir(dataPath);
ii = 0;
% rng(0);
for i = 1:length(file)
% for i = 5:8
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
			N = size(x,2);
			F = size(x,3);
			D = 2*F;	
						
			X = reshape(permute(x(1:2,:,:),[1 3 2]),D,N);
            
            % remove short trajectories
            minLen = 0;
            indShort = [];
            for kk = 1:N
                if nnz(m(kk, :)) < minLen
                    indShort = [ indShort, kk ];
                    continue;
                end
                Temp = fLabel2sLabel(m(kk, :));
                idx = (Temp(:,1)==1);
                l = Temp(idx,3) - Temp(idx,2) + 1;
                if all(l < minLen) || isempty(l)
                    indShort = [ indShort, kk ];
                end
            end
            X(:, indShort) = [];
            m(indShort, :) = [];
            s(indShort) = [];
            indX = kron(m,[1,1])';
			
% 			indX = zeros(N,D);
% 			for jj = 1:F				
% 				indX(:,2*jj-1) = m(:,jj);
% 				indX(:,2*jj) = m(:,jj);
% 			end
% 			idx_all0 = [];
%             for kk=1:N
%                 if(isequal(indX(kk,:),zeros(1,2*F)))
%                     idx_all0 = [idx_all0,kk];
%                 end
% %                 % if it is too short, remove
% %                 if nnz(indX(kk,:)) < 10
% %                     idx_all0 = [idx_all0, kk];
% %                 end
%             end
% 			indX = indX';
% 			X(:,idx_all0) =[];
% 			indX(:,idx_all0) =[];
% 			s(idx_all0) = [];			
			
			X = X.*indX;			
			
            rng('default');
            [missrate, grp, bestRank, minNcutValue] = RSIM_Incomplete(X',indX', s, 4, 1);
% 			[missrate, grp, bestRank, minNcutValue] = RSIMJBLD_Incomplete(X',indX',s,4,1);

			ii = ii+1;
			Missrate(ii) = missrate;
			disp([filePath ': ' num2str(100*Missrate(ii)) '%' ', nMotions: ' num2str(max(s)) ', seq: ' num2str(ii)]);
            clear x y s m
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
