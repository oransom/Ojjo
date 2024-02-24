function [grd_ext] = grading_extents(xyz_diff)

dat=xyz_diff;
x=dat(:,1);
y=dat(:,2);

loc = [x,y]; %create a vector set

try
    D = pdist2(loc,loc); %compute the distance between all points
catch ME
    if (strcmp(ME.identifier,'MATLAB:array:SizeLimitExceeded'))
        msg = 'Array size exceeds memory, cannot compute grading extents';
        causeException = MException('',msg);
        ME = addCause(ME,causeException);
    end
    %rethrow(ME)
    fprintf('Cannot compute grading extents - not enough memory \n');
    grd_ext=[];
    return
end

%D = pdist2(loc,loc); %compute the distance between all points
if ~isempty(D)
    try
        [idx, corepts] = dbscan(D,10,2,'Distance','precomputed'); %second variable is the search radius, confused how to automate that part
        numGroups = length(unique(idx)); %number of groups
    catch ME
        if (strcmp(ME.identifier,'MATLAB:nomem'))
            msg = 'cannot compute dbscan - out of memory';
            causeException = MException('',msg);
            ME = addCause(ME,causeException);
        end
        fprintf('cannot compute dbscan - out of memory \n');
        grd_ext=[];
        return
    end
else
    grd_ext=[];
    return
end
% figure('Visible','on')
% gscatter(loc(:,1),loc(:,2),idx,hsv(numGroups)); %plot group maps

for i=1:numGroups
    xg{i}=x(idx==i);
    yg{i}=y(idx==i);
end
for i=1:numGroups
    gbound{i}=boundary(xg{i},yg{i},1);
end
hold on
for i=1:length(gbound)
    grd_ext.gbx{i}=gbound{i};
    grd_ext.gby{i}=gbound{i};
    grd_ext.xgh{i}=xg{i};
    grd_ext.ygh{i}=yg{i};
    %plot(xgh(gbx),ygh(gby));
end
end