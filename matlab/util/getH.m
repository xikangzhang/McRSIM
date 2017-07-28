
function H = getH(features,opt)

s = size(features{1});

if ~exist('opt','var')
    opt.H_structure = 'HHt';
    opt.metric = 'JBLD';
end

H = cell(1,length(features));
for i=1:length(features)
%     t = diff(features{i},[],2);
    t = features{i};
    if strcmp(opt.H_structure,'HtH')
        Hsize = opt.H_rows;
        nc = Hsize;
        nr = size(t,1)*(size(t,2)-nc+1);
        if nr<1, error('hankel size is too large.\n'); end
        Ht = blockHankel(t,[nr nc]);
    elseif strcmp(opt.H_structure,'HHt')
        Hsize = opt.H_rows * s(1);
        nr = floor(Hsize/size(t,1))*size(t,1);
        nc = size(t,2)-floor(nr/size(t,1))+1;
        if nc<1, error('hankel size is too large.\n'); end
        Ht = blockHankel(t,[nr nc]);
    end
    H{i} = Ht;
end

end