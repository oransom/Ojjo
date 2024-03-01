function [drowh, dpileh] = ab_6_plc_nsfit_6(const,dpileh,drowh,surface)
%this function only moves the end of the effected rows, so x/y/z of the
%whole row are effected
%update all of drowh and then suck it back into drow in plc_sst;
ns_neigh_dist=const.nsdist;
%function to reduce the delta between north and south row ends


rsi=[drowh.ntpxc,drowh.ntpyc,drowh.ntpzc]; %northern points - probably can omit z
rsf=[drowh.stpxc,drowh.stpyc,drowh.stpzc]; %southern points - probably can omit z
[Idx,D]=rangesearch(rsi,rsf,ns_neigh_dist); %search for all row neighbors (n->s) within neighbor_distance
for i=1:numel(D)
    if ~isempty(D{i})
        nntr(i,1)=i; %target row (southern point of northern row)
        nnid(i,:)=Idx{i}(find(D{i}==min(D{i}))); %id of row to south (looking at its northern most point) only choose the closest one
    end
end
nntr(nntr==0)=NaN; nntr=rmmissing(nntr); %clean it all up
nnid(nnid==0)=NaN; nnid=rmmissing(nnid);

delta=drowh.stpzc(nntr)-drowh.ntpzc(nnid); %find deltas of target rows
did=find(abs(delta)>const.nsdelta); %find locations that violate
md=delta(did)-const.nsdelta; %move delta - how much the rows have to move
md=md/2; %the amount that each row should move
nntr_r=nntr(did); %nntr reduced to only violation rows;
nnid_r=nnid(did); %nnid reduced to only violation rows;
for i=1:numel(did) %add required deltas to the right places to bring everything in line
    if drowh.stpzc(nntr_r(i))>drowh.ntpzc(nnid_r(i))
        drowh.stpzc(nntr_r(i))=drowh.stpzc(nntr_r(i))-md(i);
        drowh.ntpzc(nnid_r(i))=drowh.ntpzc(nnid_r(i))+md(i);
    else
        drowh.stpzc(nntr_r(i))=drowh.stpzc(nntr_r(i))+md(i);
        drowh.ntpzc(nnid_r(i))=drowh.ntpzc(nnid_r(i))-md(i);
    end
end

id=[nntr_r;nnid_r];

%bring it back into line

np=[drowh.ntpxc(id),drowh.ntpyc(id),drowh.ntpzc(id)];
sp=[drowh.stpxc(id),drowh.stpyc(id),drowh.stpzc(id)];
uv=(sp-np)./vecnorm(sp-np,2,2); %unit vector for the slope corrected rows

for i=1:numel(id)
    dpileh.tpzc(drowh.si(id(i)))=drowh.ntpzc(id(i));
    for j=drowh.si(id(i))+1:drowh.ei(id(i)) %for each row start index+1 to end index
        span=dpileh.tpy(j)-dpileh.tpy(drowh.si(id(i))); %go back to original spans - they're in 2D
        dpileh.tpxc(j)=drowh.npx(id(i))-span.*uv(i,1); %x portion of unit vector
        dpileh.tpyc(j)=drowh.npy(id(i))-span.*uv(i,2); %y portion of unit vector
        dpileh.tpzc(j)=drowh.ntpzc(id(i))-span.*uv(i,3); %z portion of unit vector
    end
end
dpileh.bpzc=surface.F_og(dpileh.tpxc,dpileh.tpyc);

for i=1:numel(drowh.row)
    drowh.ntpxc(i)=dpileh.tpxc(drowh.si(i));
    drowh.stpxc(i)=dpileh.tpxc(drowh.ei(i));
    drowh.ntpyc(i)=dpileh.tpyc(drowh.si(i));
    drowh.stpyc(i)=dpileh.tpyc(drowh.ei(i));
    drowh.ntpzc(i)=dpileh.tpzc(drowh.si(i));
    drowh.stpzc(i)=dpileh.tpzc(drowh.ei(i));
    rise=drowh.ntpzc(i)-drowh.stpzc(i);
    run=drowh.ntpyc(i)-drowh.stpyc(i);
    drowh.slpc(i)=atand(rise./run); %if in doubt, check it out
    drowh.rowzavg(i)=mean(dpileh.tpzc(drowh.si(i):drowh.ei(i)));
end

% %% plot it to make sure
% figure
% scatter(drowh.stpxc(nntr),drowh.stpyc(nntr),'red');
% hold on
% scatter(drowh.ntpxc(nnid),drowh.ntpyc(nnid),'green')
end