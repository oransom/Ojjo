function [flp_r, flp_p] = ab_3_plc_rowflip_6(const,drowh,dpileh)
%this function moves the entire row up and down as one - therefore x's and
%y's are not effected
drowh.rowzavg_o=drowh.rowzavg;
rd=const.row_delta; %shortened version of row delta
md=const.xspacing+5;  %average row spacing in x padded by 5 feet %max distance between rows to be considered a neighbor
%%initial loop
cbf_o=[]; %set to nothing so it'll run once
oos=[];
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
for i=1:numel(oos) %doing this in a loop for clarity - could be vectorized - flip additional rows due to exposure
    if ~isnan(drowh.nne(oos(i)))
        zdel=drowh.rowzavg(oos(i))-drowh.rowzavg(drowh.nne(oos(i)));
        if zdel < -rd %if row to the east of an oos row is causing it to flip by being too high it becomes an exterior
            flp_r.flip2ext(drowh.nne(oos(i)))=1;
            if ~isnan(drowh.nne(drowh.nne(oos(i))))
                zdel=drowh.rowzavg(oos(i))-drowh.rowzavg(drowh.nne(drowh.nne(oos(i))));
                if zdel < -rd %if the row to the row to the east of an oos row is causing it to flip by being too high it becomes an exterior
                    flp_r.flip2ext(drowh.nne(drowh.nne(oos(i))))=1;
                end
            end
        end
    end
    if ~isnan(drowh.nnw(oos(i)))
        zdel=drowh.rowzavg(oos(i))-drowh.rowzavg(drowh.nnw(oos(i)));
        if zdel < -rd %if row to the east of an oos row is causing it to flip by being too high it becomes an exterior
            flp_r.flip2ext(drowh.nnw(oos(i)))=1;
            if ~isnan(drowh.nnw(drowh.nnw(oos(i))))
                zdel=drowh.rowzavg(oos(i))-drowh.rowzavg(drowh.nnw(drowh.nnw(oos(i))));
                if zdel < -rd %if the row to the row to the east of an oos row is causing it to flip by being too high it becomes an exterior
                    flp_r.flip2ext(drowh.nnw(drowh.nnw(oos(i))))=1;
                end
            end
        end
    end
end
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
