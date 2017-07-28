function velocity = getVelocity(features, step)
% get features' first order derivative along the time dimension
if nargin < 2
    step = 1;
end
h = [1 zeros(1,step-1) -1];

velocity = cell(1,length(features));
for i=1:length(features)
    velocity{i} = conv2(features{i}, h,'valid');
end

end