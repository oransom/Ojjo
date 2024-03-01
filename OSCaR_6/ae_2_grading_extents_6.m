function [surface] = ae_2_grading_extents_6(surface)


xh=reshape(surface.xq,[],1);
yh=reshape(surface.yq,[],1);
zh=reshape(surface.sdiff,[],1);
x=xh(~isnan(zh));
y=yh(~isnan(zh));

loc = [x,y]; %create a vector set

if isempty(loc) | height(loc)>10000 %no/don't grade
    surface.xg=[];
    surface.yg=[];
    surface.gbound=[];
    return
end

% try
%     D = pdist2(loc,loc); %compute the distance between all points
% catch ME
%     if (strcmp(ME.identifier,'MATLAB:array:SizeLimitExceeded'))
%         msg = 'Array size exceeds memory, cannot compute grading extents';
%         causeException = MException('',msg);
%         ME = addCause(ME,causeException);
%     end
%     %rethrow(ME)
%     fprintf('Cannot compute grading extents - not enough memory \n');
%     surface.gbound=[];
%     return
% end

%D = pdist2(loc,loc); %compute the distance between all points
%if ~isempty(D)
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
%hold on
% for i=1:length(gbound)
%     surface.gbound=gbound;
%     surface.grd_ext_gby=gbound{i};
%     surface.grd_ext_gbz=gbound{i};
%     surface.grd_ext_xgh{i,1}=xg{i};
%     surface.grd_ext_ygh{i,1}=yg{i};
%     surface.grd_ext_zgh{i,1}=ones(numel(yg{i}),1).*100000;
%     %plot(xgh(gbx),ygh(gby));
% end
end