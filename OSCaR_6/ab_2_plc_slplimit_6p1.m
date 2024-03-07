function [drowh, dpileh] = ab_2_plc_slplimit_6p1(const,drowh,dpileh)
%function to limit north and south row slope to allowable levels
%find first and last tpzf in each row
%north facing slope is negative
%keeping north and south separated for readability
idn=find(drowh.slp<0 & drowh.slp<-const.nslopelimit); %find rows exceeding north slope limit
ids=find(drowh.slp>0 & drowh.slp>const.sslopelimit); %find rows exceeding south slope limit
id=sort([idn;ids]);
ns_diff=drowh.slp(idn)+const.nslopelimit; %find delta between limit and measured slope
ss_diff=drowh.slp(ids)-const.sslopelimit;
ns_ydiff=dpileh.tpyc(drowh.si(idn))-dpileh.tpyc(drowh.ei(idn)); %row length at this stage
ss_ydiff=dpileh.tpyc(drowh.si(ids))-dpileh.tpyc(drowh.ei(ids)); %row length at this stage
zn_diff_calc=tand(ns_diff).*ns_ydiff; %calculate vertical differential
zs_diff_calc=tand(ss_diff).*ss_ydiff;

slp.tpx=dpileh.tpxc; %preallocating holder array
slp.tpy=dpileh.tpyc;
slp.tpz=dpileh.tpzc;

for i=1:numel(idn)
    slp.tpz(drowh.si(idn(i)))=slp.tpz(drowh.si(idn(i)))+zn_diff_calc(i)/2; %make rows face more south
    slp.tpz(drowh.ei(idn(i)))=slp.tpz(drowh.ei(idn(i)))-zn_diff_calc(i)/2;
end
for i=1:numel(ids)
    slp.tpz(drowh.si(ids(i)))=slp.tpz(drowh.si(ids(i)))-zs_diff_calc(i)/2; %make rows face more north
    slp.tpz(drowh.ei(ids(i)))=slp.tpz(drowh.ei(ids(i)))+zs_diff_calc(i)/2;
end
%check
rise=slp.tpz(drowh.si)-slp.tpz(drowh.ei);
run=slp.tpy(drowh.si)-slp.tpy(drowh.ei);
slp.scheck=atand(rise./run); %if in doubt, check it out

%find unit vector from first to last pile of adjusted rows and realign
%whole row
np=[slp.tpx(drowh.si(id)),slp.tpy(drowh.si(id)),slp.tpz(drowh.si(id))];
sp=[slp.tpx(drowh.ei(id)),slp.tpy(drowh.ei(id)),slp.tpz(drowh.ei(id))];
uv=(sp-np)./vecnorm(sp-np,2,2); %unit vector for the slope corrected rows

for i=1:numel(id)
    for j=drowh.si(id(i))+1:drowh.ei(id(i)) %for each row start index+1 to end index
        span=dpileh.tpy(j)-dpileh.tpy(drowh.si(id(i))); %go back to original spans - they're in 2D
        slp.tpx(j)=drowh.npx(id(i))-span.*uv(i,1); %x portion of unit vector
        slp.tpy(j)=drowh.npy(id(i))-span.*uv(i,2); %y portion of unit vector
        slp.tpz(j)=drowh.npz(id(i))-span.*uv(i,3); %z portion of unit vector
    end
end

dpileh.tpxc=slp.tpx;
dpileh.tpyc=slp.tpy;
dpileh.tpzc=slp.tpz;
drowh.slpc=slp.scheck;
drowh.rlengthc=drowh.npy-drowh.spyc;
drowh.ntpxc=dpileh.tpxc(drowh.si); %northern tpx corrected
drowh.stpxc=dpileh.tpxc(drowh.ei); %southern tpx corrected
drowh.ntpyc=dpileh.tpyc(drowh.si); %northern tpy corrected
drowh.stpyc=dpileh.tpyc(drowh.ei); %southern tpy corrected
drowh.ntpzc=dpileh.tpzc(drowh.si); %northern selected tpz corrected
drowh.stpzc=dpileh.tpzc(drowh.ei); %southern selected tpz corrected

end