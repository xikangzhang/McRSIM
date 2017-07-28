function VR = findVR(X, mask, UpperR)

[N,F] = size(X);

Ns = sum(mask);
[Ns_sort,idx_sort] = sort(Ns,'Descend');
idx_init = idx_sort(:,1:UpperR);
mask_init = mask(:,idx_init);
sum_mask_init = sum(mask_init,2);
idx_delete_row = find(sum_mask_init<UpperR);

Xinit = X(:,idx_init);

Xinit(idx_delete_row,:) = []; % delete these row with missing data

for i=1:N-length(idx_delete_row)
	if(isequal(Xinit(i,:),zeros(1,UpperR)))
		Xinit(i,:) = Xinit(i-1,:);
	end
end

[Vinit,~,~] = svd(Xinit,'econ');
Vinit = Vinit(:,1:UpperR);

V = zeros(N,UpperR);
for i=1:N
	if(isempty(find(idx_delete_row==i)))
		aug_idx = sort([i;idx_delete_row]);
		V(i,:) = Vinit(i-find(aug_idx==i)+1,:);
	end
end
Vinit = V;

idx = find(mask(:)==1);
[I,J] = ind2sub([N,F], idx);
X = X(idx);

step_size = .01;
maxCycles = 200;

[VR] = grouse(I,J,X,N,F,UpperR,step_size,maxCycles,Vinit);



end