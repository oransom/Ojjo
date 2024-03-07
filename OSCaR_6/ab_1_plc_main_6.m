function [drow, dpile]=ab_1_plc_main_6(const, surface, dpile, drow)

switch char(const.rsolve)
    case 'slope' %slope solution
        for i=1:height(drow)
            rhld.pi=[drow.npx,drow.npy,drow.npz+const.min_wp]; %pile initial
            rhld.pf=[drow.npx,drow.spyc,drow.spz+const.min_wp]; %pile final
            rhld.npz=drow.npz+const.min_wp; %holder for npz
        end
    case 'rbestfit' %row best fit solution
        for i=1:height(drow)
            rhld.pi=[drow.npx,drow.npy,drow.lbf_npz+const.min_wp];
            rhld.pf=[drow.npx,drow.spyc,drow.lbf_spz+const.min_wp];
            rhld.npz=drow.lbf_npz+const.min_wp;
        end
    case 'nbestfit' %neighbor best fit solution
        for i=1:height(drow)
            rhld.pi=[drow.npx,drow.npy,drow.nbf_npz+const.min_wp];
            rhld.pf=[drow.npx,drow.spyc,drow.nbf_spz+const.min_wp];
            rhld.npz=drow.nbf_npz+const.min_wp;
        end
end
drow = removevars(drow,{'lbf_npz','lbf_spz','nbf_npz','nbf_spz'});
%find unit vector from first to last pile
uv=(rhld.pf-rhld.pi)./vecnorm(rhld.pf-rhld.pi,2,2); %unit vector for any of the straight tt sol's above
dpile.FFR=zeros(numel(dpile.tpx),1); %FLAG FOR REVIEW PILE
dpile.tpxc=dpile.tpx; dpile.tpyc=dpile.tpy; dpile.tpxc=dpile.tpz;
for i=1:height(drow)
    dpile.tpxc(drow.si(i))=drow.npx(i); %set the first pile in each row
    dpile.tpyc(drow.si(i))=drow.npy(i);
    dpile.tpzc(drow.si(i))=rhld.npz(i);
    for j=drow.si(i)+1:drow.ei(i) %for each row start index+1 to end index
        span=dpile.tpy(j)-dpile.tpy(drow.si(i));
        dpile.tpxc(j)=drow.npx(i)-span.*uv(i,1); %x portion of unit vector
        dpile.tpyc(j)=drow.npy(i)-span.*uv(i,2); %y portion of unit vector
        dpile.tpzc(j)=rhld.npz(i)-span.*uv(i,3); %z portion of unit vector
    end
end
clear span uv rhld
for i=1:height(drow) %calculate slope - negative slope faces north
    drow.slp(i)=atand((dpile.tpzc(drow.ei(i))-dpile.tpzc(drow.si(i)))/...
        (dpile.tpyc(drow.ei(i))-dpile.tpyc(drow.si(i))));
end
%slope limiting - done before minimum pile heights or flood pile heights
%this is done agnostic of the solution type above
%dpile.tpx/y/z c values are overwritten
%drow.slp is overwritten, but np/sp values in drow are not - will be
%cleaned up after this function
if const.slopelimit==1
    [drow, dpile] = ab_2_plc_slplimit_6p1(const,drow,dpile);
else
    drow.slpc=drow.slp;
end

switch char(const.rsolve)
    case 'slope' %slope solution
        phld.ih=dpile.tpzc-surface.F_og(dpile.tpxc,dpile.tpyc);
        for i=1:height(drow) %ensure that no pile is less than minimum required pile
            j=drow.si(i):drow.ei(i);
            rhld.mh(i,1)=min(phld.ih(j));
            rhld.add(i,1)=const.min_wp-rhld.mh(i);
            dpile.tpzc(j)=dpile.tpzc(j)+rhld.add(i); %tpzc2 is post-fit solution
        end
        clear phld rhld
    case 'rbestfit' %row best fit solution
        %do not change the initial best fit solution
    case 'nbestfit' %neighbor best fit solution
        %do not change the initial best fit solution
end

%ns align should occur here, after rows are moved above
if const.nsfit==1
    [drow, dpile] = ab_6_plc_nsfit_6p1(const,dpile,drow,surface);
end

%% Flood calcs
if ~strcmpi(const.flood,'na') %flood surface should always be above ground surface
    fld2grd=surface.Flood_s(dpile.tpxc,dpile.tpyc)-surface.F_og(dpile.tpxc,dpile.tpyc);
    fld2grd(fld2grd<0.01)=0;
    cottd=const.min_wp-const.freeboard; %if wp is above freeboard it'll be positive
    phld.fh=dpile.tpzc-surface.Flood_s(dpile.tpxc,dpile.tpyc); %top of pile z to flood surface 
    for i=1:height(drow) %ensure that no pile is less than minimum required pile
        j=drow.si(i):drow.ei(i);
        tpzcf=surface.Flood_s(dpile.tpxc(drow.si(i)-1+find(fld2grd(j)>0)),dpile.tpyc(drow.si(i)-1+find(fld2grd(j)>0)))+const.freeboard; %find wet top of pile in row
        tpzd=tpzcf-dpile.tpzc(drow.si(i)-1+find(fld2grd(j)>0));
        %minfchk(i,1)=const.freeboard-min(dpile.tpzc(j)-surface.Flood_s(dpile.tpxc(j),dpile.tpyc(j)));%trying to sort out an error
        if max(tpzd)>0
            rhld.fadd(i,1)=max(tpzd(tpzd>=0));
        else
            rhld.fadd(i,1)=0;     
        end
        dpile.tpzc(j)=dpile.tpzc(j)+rhld.fadd(i);
        dpile.matpzc(j)=dpile.tpzc(j)-surface.F_og(dpile.tpxc(j),dpile.tpyc(j)); %minimum allowable tpzc
    end
    drow.ntpxc=dpile.tpxc(drow.si); %northern tpx corrected
    drow.stpxc=dpile.tpxc(drow.ei); %southern tpx corrected
    drow.ntpyc=dpile.tpyc(drow.si); %northern tpy corrected
    drow.stpyc=dpile.tpyc(drow.ei); %southern tpy corrected
    drow.ntpzc=dpile.tpzc(drow.si); %northern selected tpz corrected
    drow.stpzc=dpile.tpzc(drow.ei); %southern selected tpz corrected
else
    %tpzc is unchanged since no flood
    drow.ntpxc=dpile.tpxc(drow.si); %northern tpx corrected
    drow.stpxc=dpile.tpxc(drow.ei); %southern tpx corrected
    drow.ntpyc=dpile.tpyc(drow.si); %northern tpy corrected
    drow.stpyc=dpile.tpyc(drow.ei); %southern tpy corrected
    drow.ntpzc=dpile.tpzc(drow.si); %northern selected tpz corrected
    drow.stpzc=dpile.tpzc(drow.ei); %southern selected tpz corrected
end

dpile.bpzc=surface.F_og(dpile.tpxc,dpile.tpyc);
for i=1:height(drow) %calculate how much you can move each row up - for flipex calcs later
    j=drow.si(i):drow.ei(i);
    drow.prmng(i)=const.max_wp-max(dpile.tpzc(j)-dpile.bpzc(j));
end
clear phld rhld i j

for i=1:height(drow)
    drow.rowzavg(i)=mean(dpile.tpzc(drow.si(i):drow.ei(i)));
end
clear i

if const.flipex==1 %flip rows to exterior based on ATI and NXT best doc 
    switch char(const.tracker)
        case {'ATI','Ojjo_ATI'} %ATI Row Flip
            [flp_r, flp_p] = ab_3_plc_rowflip_6(const,drow,dpile); %returns flip row and flip pile data
            %roll back into drow and dpile
            drow.rowzavg=flp_r.nrowzavg;
            drow.flip2ext=flp_r.flip2ext;
            drow.ntpzc=drow.ntpzc+flp_r.add2row;
            drow.stpzc=drow.stpzc+flp_r.add2row;
            dpile.tpzc=flp_p.tpz;
            clear flp_p flp_r
        case {'NXT','Ojjo_NXT'} %NXT Row Flip
            [flp_red, flp_ped] = ab_4_plc_rowflip_6_nxt_edge(const,drow,dpile); %returns flip row and flip pile data
            %roll back into drow and dpile
            drow.rowzavg=flp_red.nrowzavg;
            drow.flip2edg=flp_red.flip2ext;
            drow.ntpzc=drow.ntpzc+flp_red.add2row;
            drow.stpzc=drow.stpzc+flp_red.add2row;
            dpile.tpzc=flp_ped.tpz;
            clear flp_ped flp_red
            [flp_rex, flp_pex] = ab_5_plc_rowflip_6_nxt_ext(const,drow,dpile); %returns flip row and flip pile data
            %roll back into drow and dpile
            drow.rowzavg=flp_rex.nrowzavg;
            drow.flip2ext=flp_rex.flip2ext;
            drow.ntpzc=drow.ntpzc+flp_rex.add2row;
            drow.stpzc=drow.stpzc+flp_rex.add2row;
            dpile.tpzc=flp_pex.tpz;
            clear flp_pex flp_rex
            drow.flip2edg(drow.flip2ext==1)=0;
    end
end

%%clean and arrange table
drow = removevars(drow,{'spyc','spzc'});
drow = movevars(drow, "block", "Before", "npx");
drow = movevars(drow, "sect", "Before", "npx");
drow = movevars(drow, "row", "Before", "npx");
drow = movevars(drow, "motorloc", "Before", "npx");
drow = movevars(drow, "si", "Before", "motorloc");
drow = movevars(drow, "ei", "Before", "motorloc");
drow = movevars(drow, "span", "Before", "npx");
drow = movevars(drow, "npz", "Before", "mpy");
drow = movevars(drow, "mdptx", "Before", "mpy");
drow = movevars(drow, "mdpty", "Before", "mpy");
drow = movevars(drow, "spz", "Before", "nnw");

%align motors
if const.motorc==1 %tpyc2 is motor adjusted - if selected, assume z doesn't shift - maybe an issue?
    for i=1:height(drow)
        j=drow.si(i):drow.ei(i);
        m=drow.si(i)-1+drow.motorloc(i); %location of motor in row in whole dpile table
        y_adj=dpile.tpyc(m)-dpile.tpy(m); %adjustment needed between non-slope and solved
        dpile.tpyc(j)=dpile.tpyc(j)-y_adj; %adjust row back to its original motor location
    end
end

if contains(char(const.tracker),'ojjo','IgnoreCase',true)
    [drow, dpile]=ab_7_plc_ojjo_6(const, surface, dpile, drow);
elseif strcmpi(char(const.tracker),'ATI') | strcmpi(char(const.tracker),'NXT')
    [drow, dpile]=ab_8_plc_atinxt_6(const, surface, dpile, drow);
end

%% check with plotting
% trow=1;
% er_diff=dpile.tpy(drow.ei(trow))-dpile.tpyc(drow.ei(trow));
% figure
% scatter(dpile.tpy(drow.si(trow):drow.ei(trow)),dpile.tpzc(drow.si(trow):drow.ei(trow)),'red');
% hold on
% plot(dpile.tpyc(drow.si(trow):drow.ei(trow)),dpile.tpzc(drow.si(trow):drow.ei(trow)))
% scatter(dpile.tpyc(drow.si(trow):drow.ei(trow)),dpile.tpzc(drow.si(trow):drow.ei(trow)),'green');

end


