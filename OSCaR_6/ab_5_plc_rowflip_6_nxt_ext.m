function [flp_r, flp_p] = ab_5_plc_rowflip_6_nxt_ext(const,drowh,dpileh)
drowh.rowzavg_o=drowh.rowzavg;
int2ext=0.12;
md=const.xspacing+5; 
rd=int2ext*const.xspacing; %row delta is calced as the mode of xrow spacing
%%initial loop
oos=[];
cbf_o=[]; %set to nothing so it'll run once
[cbf, ~, ~, ~, drowh] = exflip(drowh, rd, md); %get initials back for can be fixed
while ~isequal(cbf,cbf_o) %run until no more can be fixed
    [cbf, did, del_r, ~, drowh] = exflip(drowh, rd, md);
    drowh.rowzavg(cbf)=drowh.rowzavg(cbf)+(del_r(did)-rd); %increase the rowzavg of the rows that can be moved
    drowh.prmng(cbf)=drowh.prmng(cbf)-(del_r(did)-rd); %decrease the pile remaining of the piles that were moved
    [cbf_o, ~, ~, oos, ~] = exflip(drowh, rd, md); %checker for cbf_o and final on oos
end
%%flip2ex required adjacent rows
flp_r.nrowzavg=drowh.rowzavg; %return of new rowzavg
flp_r.flip2ext=zeros(numel(drowh.row),1);
flp_r.flip2ext(oos)=1; %flip to exteriors
flp_r.add2row=drowh.rowzavg-drowh.rowzavg_o;
drowh.rowzavg=flp_r.nrowzavg;
drowh.flip2ext=flp_r.flip2ext;
drowh.ntpzc=drowh.ntpzc+flp_r.add2row;
drowh.stpzc=drowh.stpzc+flp_r.add2row;
id = find(flp_r.add2row);
flp_p.tpz=dpileh.tpzc;

for i=1:numel(id)
j=drowh.si(id(i))+1:drowh.ei(id(i)); %for each row start index+1 to end index
flp_p.tpz(j)=flp_p.tpz(j)+flp_r.add2row(i); %z is moved uniformly up
end
end

function [cbf, did, del_r, oos, drowh] = exflip(drowh, rd, md)
wd=[]; wv=[]; ed=[]; ev=[]; %so no errors if they don't exist
for i=1:height(drowh)
    if ~isnan(drowh.nnw(i))
        if drowh.rowzavg(drowh.nnw(i))-drowh.rowzavg(i)>rd &&...
                abs(drowh.ntpxc(drowh.nnw(i))-drowh.ntpxc(i))<md%west is too tall and close enough
            wd(i,1)=drowh.rowzavg(drowh.nnw(i))-drowh.rowzavg(i);
            wv(i,1)=i;
        end
    end
    if ~isnan(drowh.nne(i))
        if drowh.rowzavg(drowh.nne(i))-drowh.rowzavg(i)>rd &&...
                abs(drowh.ntpxc(drowh.nne(i))-drowh.ntpxc(i))<md %east is too tall and close enough
            ed(i,1)=drowh.rowzavg(drowh.nne(i))-drowh.rowzavg(i);
            ev(i,1)=i;
        end
    end
end

oos=[wv;ev]; %put them in a single array
oos(oos==0)=NaN;
oos=rmmissing(oos); %strip array to only oos
del_r=[wd;ed]; %put them in a single array
del_r(del_r==0)=NaN;
del_r=rmmissing(del_r); %strip array to only oos
cbf=oos((del_r-rd)<drowh.prmng(oos)); %cbf = can be fixed
did=find((del_r-rd)<drowh.prmng(oos)); %locations in del_r array of cbf
end
