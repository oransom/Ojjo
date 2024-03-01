function [slp, id] = ab_2_plc_slplimit_6(const,drow,dpile)
%function to limit north and south row slope to allowable levels
%find first and last tpzf in each row
%north facing slope is negative
%keeping north and south separated for readability
idn=find(drow.slp<0 & drow.slp<-const.nslopelimit); %find rows exceeding north slope limit
ids=find(drow.slp>0 & drow.slp>const.sslopelimit); %find rows exceeding south slope limit
id=sort([idn;ids]);
ns_diff=drow.slp(idn)+const.nslopelimit; %find delta between limit and measured slope
ss_diff=drow.slp(ids)-const.sslopelimit;
ns_ydiff=dpile.tpyc(drow.si(idn))-dpile.tpyc(drow.ei(idn)); %row length at this stage
ss_ydiff=dpile.tpyc(drow.si(ids))-dpile.tpyc(drow.ei(ids)); %row length at this stage
zn_diff_calc=tand(ns_diff).*ns_ydiff; %calculate vertical differential
zs_diff_calc=tand(ss_diff).*ss_ydiff;

slp.tpx=dpile.tpxc; %preallocating holder array
slp.tpy=dpile.tpyc;
slp.tpz=dpile.tpzc;

for i=1:numel(idn)
    slp.tpz(drow.si(idn(i)))=slp.tpz(drow.si(idn(i)))+zn_diff_calc(i)/2; %make rows face more south
    slp.tpz(drow.ei(idn(i)))=slp.tpz(drow.ei(idn(i)))-zn_diff_calc(i)/2;
end
for i=1:numel(ids)
    slp.tpz(drow.si(ids(i)))=slp.tpz(drow.si(ids(i)))-zs_diff_calc(i)/2; %make rows face more north
    slp.tpz(drow.ei(ids(i)))=slp.tpz(drow.ei(ids(i)))+zs_diff_calc(i)/2;
end
%check
rise=slp.tpz(drow.si(id))-slp.tpz(drow.ei(id));
run=slp.tpy(drow.si(id))-slp.tpy(drow.ei(id));
slp.scheck=atand(rise./run); %if in doubt, check it out

%find unit vector from first to last pile of adjusted rows and realign
%whole row
np=[slp.tpx(drow.si(id)),slp.tpy(drow.si(id)),slp.tpz(drow.si(id))];
sp=[slp.tpx(drow.ei(id)),slp.tpy(drow.ei(id)),slp.tpz(drow.ei(id))];
uv=(sp-np)./vecnorm(sp-np,2,2); %unit vector for the slope corrected rows

for i=1:numel(id)
    for j=drow.si(id(i))+1:drow.ei(id(i)) %for each row start index+1 to end index
        span=dpile.tpy(j)-dpile.tpy(drow.si(id(i))); %go back to original spans - they're in 2D
        slp.tpx(j)=drow.npx(id(i))-span.*uv(i,1); %x portion of unit vector
        slp.tpy(j)=drow.npy(id(i))-span.*uv(i,2); %y portion of unit vector
        slp.tpz(j)=drow.npz(id(i))-span.*uv(i,3); %z portion of unit vector
    end
end
end