
function [HH, H] = getHH_missing(X, mask, opt)

if ~exist('opt','var')
    opt.H_structure = 'HHt';
    opt.metric = 'JBLD';
end

n = size(X, 2);
if nargout > 1
    H = cell(1,n);
end

HH = cell(1,n);
for i = 1:n
    t = reshape(X(:, i), 2, []);
    ind = reshape(mask(:, i), 2, []);
    ind = ind(1, :);
    seg = trajSplit(t, ind);
    seg2 = removeShortTraj(seg, opt.H_rows);
    Hs = getH(seg2, opt);
    if strcmp(opt.H_structure, 'HHt')
        Ht = cell2mat(Hs);
        HHt = Ht * Ht';
    elseif strcmp(opt.H_structure, 'HtH')
        Ht = cell2mat(Hs');
        HHt = Ht' * Ht;
    end
    
    HHt = HHt / norm(HHt, 'fro');
    if strcmp(opt.metric,'JBLD') ...
            || strcmp(opt.metric,'AIRM') || strcmp(opt.metric,'LERM')...
            || strcmp(opt.metric,'KLDM')
        I = opt.sigma*eye(size(HHt));
        HH{i} = HHt + I;
    elseif strcmp(opt.metric,'binlong') || strcmp(opt.metric,'SubspaceAngle') ||...
            strcmp(opt.metric,'SubspaceAngleFast')
        HH{i} = HHt;
    end
    
    if nargout > 1
        H{i} = Ht;
    end
end

end