function [Missrate, index] = ErrorRate2(Seg,RefSeg)

if(max(max(RefSeg))~=1 && max(max(Seg))~=1) %Calculate the error rate using the results of standard spectral clustering method
	%Missrate = Misclassification(Seg,RefSeg);
    [miss, index] = missclassGroups( Seg,RefSeg,max(RefSeg) );
	Missrate = miss ./ length(RefSeg);
elseif(max(max(RefSeg))==1 && max(max(Seg))==1)
	[N,n] = size(RefSeg); %number of clusters
	SegVec = zeros(N,1);
	RefSegVec = zeros(N,1);
	for i=1:n
		SegVec = SegVec+i*Seg(:,i);
		RefSegVec = RefSegVec+i*RefSeg(:,i);
	end
	%Missrate = Misclassification(SegVec,RefSegVec);
    [miss, index] = missclassGroups( SegVec,RefSegVec,max(RefSegVec) );
	Missrate = miss ./ length(RefSegVec);
elseif(max(max(RefSeg))==1 && max(max(Seg))~=1)
	[N,n] = size(RefSeg); %number of clusters
	RefSegVec = zeros(N,1);
	for i=1:n
		RefSegVec = RefSegVec+i*RefSeg(:,i);		
	end
	%Missrate = Misclassification(Seg,RefSegVec);
    [miss, index] = missclassGroups( Seg,RefSegVec,max(RefSegVec) );
	Missrate = miss ./ length(RefSegVec);
else
	[N,n] = size(Seg); %number of clusters
	SegVec = zeros(N,1);
	for i=1:n
		SegVec = SegVec+i*Seg(:,i);		
	end
	%Missrate = Misclassification(SegVec,RefSeg);
    [miss, index] = missclassGroups( SegVec,RefSeg,max(RefSeg) );
	Missrate = miss ./ length(RefSeg);
end
