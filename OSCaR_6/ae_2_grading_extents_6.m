function [surface] = ae_2_grading_extents_6(surface)


xh=reshape(surface.xq,[],1);
yh=reshape(surface.yq,[],1);
zh=reshape(surface.sdiff,[],1);
x=xh(~isnan(zh));
y=yh(~isnan(zh));

loc = [x,y]; %create a vector set

if isempty(loc) | height(loc)>1000000 %no/don't grade
    surface.xg=[];
    surface.yg=[];
    surface.gbound=[];
    return
end


    try
        [idx, ~] = dbscan(loc,25,2);%,'Distance','precomputed'); %second variable is the search radius, confused how to automate that part
        numGroups = length(unique(idx)); %number of groups
    catch ME
        if (strcmp(ME.identifier,'MATLAB:nomem'))
            msg = 'cannot compute dbscan - out of memory';
            causeException = MException('',msg);
            ME = addCause(ME,causeException);
        end
        fprintf('cannot compute dbscan - out of memory \n');
        surface.gbound=[];
        return
    end
% else
%     grd_ext=[];
%     return
% end
% figure('Visible','on')
% gscatter(loc(:,1),loc(:,2),idx,hsv(numGroups)); %plot group maps

for i=1:numGroups
    surface.xg{i}=x(idx==i);
    surface.yg{i}=y(idx==i);
    surface.gbound{i}=boundary(surface.xg{i},surface.yg{i},1);
end
end