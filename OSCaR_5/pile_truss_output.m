function [po, to, pmd, lengthst] = pile_truss_output(kpv, const, cpl, trs, pmd)

for i=1:length(kpv)
    po.(kpv{i}).tpx=reshape(cpl.(kpv{i}).tpx,[],1);
    po.(kpv{i}).tpy=reshape(cpl.(kpv{i}).tpy,[],1);
    po.(kpv{i}).tpz=reshape(cpl.(kpv{i}).tpzf,[],1);
    po.(kpv{i}).t=[po.(kpv{i}).tpx,po.(kpv{i}).tpy,po.(kpv{i}).tpz];
    po.(kpv{i}).bpx=reshape(cpl.(kpv{i}).bpx,[],1);
    po.(kpv{i}).bpy=reshape(cpl.(kpv{i}).bpy,[],1);
    po.(kpv{i}).bpz=reshape(cpl.(kpv{i}).bpzf,[],1);
    po.(kpv{i}).nsy=reshape(cpl.(kpv{i}).nsy,[],1);
    po.(kpv{i}).nsx=reshape(cpl.(kpv{i}).nsx,[],1);
    po.(kpv{i}).b=[po.(kpv{i}).bpx,po.(kpv{i}).bpy,po.(kpv{i}).bpz];
    po.(kpv{i}).slp=reshape(cpl.(kpv{i}).slp',[],1);
    po.(kpv{i}).type=reshape(cpl.(kpv{i}).piletype,[],1);
    if const.Nevados==1
        po.(kpv{i}).upslope=reshape(cpl.(kpv{i}).slopeu',[],1);
        po.(kpv{i}).downslope=reshape(cpl.(kpv{i}).sloped',[],1);
        po.(kpv{i}).pierSumSlope=reshape(cpl.(kpv{i}).sumslope',[],1);
        po.(kpv{i}).NevJoint=reshape(cpl.(kpv{i}).NevJoint',[],1);
    end

    to.(kpv{i}).wx=reshape(trs.(kpv{i}).wx,[],1);
    to.(kpv{i}).wy=reshape(trs.(kpv{i}).wy,[],1);
    to.(kpv{i}).wz=reshape(trs.(kpv{i}).wz,[],1);
    to.(kpv{i}).w=[to.(kpv{i}).wx,to.(kpv{i}).wy,to.(kpv{i}).wz];
    to.(kpv{i}).ex=reshape(trs.(kpv{i}).ex,[],1);
    to.(kpv{i}).ey=reshape(trs.(kpv{i}).ey,[],1);
    to.(kpv{i}).ez=reshape(trs.(kpv{i}).ez,[],1);
    to.(kpv{i}).e=[to.(kpv{i}).ex,to.(kpv{i}).ey,to.(kpv{i}).ez];
    if const.thirdleg==1
        to.(kpv{i}).mx=reshape(trs.(kpv{i}).mx,[],1);
        to.(kpv{i}).my=reshape(trs.(kpv{i}).my,[],1);
        to.(kpv{i}).mz=reshape(trs.(kpv{i}).mz,[],1);
        to.(kpv{i}).m=[to.(kpv{i}).mx,to.(kpv{i}).my,to.(kpv{i}).mz];
        to.(kpv{i}).lm=zeros(length(po.(kpv{i}).t),1); %have to do this to ensure same length as le and lw
    end
    for j=1:length(po.(kpv{i}).t)
        po.(kpv{i}).lp(j)=pdist2(po.(kpv{i}).t(j,:),po.(kpv{i}).b(j,:));
        to.(kpv{i}).le(j)=pdist2(po.(kpv{i}).t(j,:),to.(kpv{i}).e(j,:));
        to.(kpv{i}).lw(j)=pdist2(po.(kpv{i}).t(j,:),to.(kpv{i}).w(j,:));
        if const.thirdleg==1
            if to.(kpv{i}).m(j) ~= 0
                to.(kpv{i}).lm(j)=pdist2(po.(kpv{i}).t(j,:),to.(kpv{i}).m(j,:));
            end
        end
    end
    to.(kpv{i}).le=reshape(to.(kpv{i}).le,[],1);
    to.(kpv{i}).lw=reshape(to.(kpv{i}).lw,[],1);
    if const.thirdleg==1
        to.(kpv{i}).lm=reshape(to.(kpv{i}).lm,[],1);
    end
    po.(kpv{i}).lp=reshape(po.(kpv{i}).lp,[],1);
    lengths.(kpv{i}).out=[to.(kpv{i}).lw,po.(kpv{i}).lp,to.(kpv{i}).le];
    po.(kpv{i}).pot=struct2table(po.(kpv{i}));
    lengthst.(kpv{i}).out=struct2table(lengths.(kpv{i}));

    k=0;
    for j=1:length(po.(kpv{i}).bpx)
        if ~(isnan(po.(kpv{i}).bpx(j)))
            k=k+1;
            pmd.(kpv{i}).fn(j)=k;
        else
            pmd.(kpv{i}).fn(j)=0;
        end
    end
    pmd.(kpv{i}).fn=pmd.(kpv{i}).fn';
    for j=1:length(po.(kpv{i}).bpx)
        pmd.(kpv{i}).fn_s{j}=[kpv{i} '_' num2str(pmd.(kpv{i}).fn(j))];
        pmd.(kpv{i}).block{j}=[kpv{i}];
    end
    for j=1:length(po.(kpv{i}).bpx)
        pmd.(kpv{i}).id_number{j}=[pmd.(kpv{i}).tl{j} '' num2str(pmd.(kpv{i}).sect(j)) '' num2str(pmd.(kpv{i}).row(j)*100)];
        pmd.(kpv{i}).bsr{j}=[(kpv{i}) '.' num2str(pmd.(kpv{i}).sect(j)) '.' num2str(pmd.(kpv{i}).row(j))]; %make block section row identifier for slope issue
        pmd.(kpv{i}).bsrt{j}=[(kpv{i}) '.' num2str(pmd.(kpv{i}).sect(j)) '.' num2str(pmd.(kpv{i}).row(j)) '.' pmd.(kpv{i}).tl{j}];
    end
    po.(kpv{i}).lpmax=zeros(length(po.(kpv{i}).lp),1);
    po.(kpv{i}).lpmin=zeros(length(po.(kpv{i}).lp),1);
    po.(kpv{i}).lpin=zeros(length(po.(kpv{i}).lp),1);

    for j=1:length(po.(kpv{i}).lp)
        if po.(kpv{i}).lp(j)<cpl.(kpv{i}).minlenf(j)-0.02
            po.(kpv{i}).lpmin(j)=po.(kpv{i}).lp(j);
        elseif po.(kpv{i}).lp(j)>const.max_wp+0.02
            po.(kpv{i}).lpmax(j)=po.(kpv{i}).lp(j);
        else
            po.(kpv{i}).lpin(j)=po.(kpv{i}).lp(j);
        end
    end

    po.(kpv{i}).lpmin(po.(kpv{i}).lpmin==0)=NaN;
    po.(kpv{i}).lpmax(po.(kpv{i}).lpmax==0)=NaN;
    po.(kpv{i}).lpin(po.(kpv{i}).lpin==0)=NaN;
    po.(kpv{i}).slpln=zeros(1,length(po.(kpv{i}).slp));
    po.(kpv{i}).slpn=zeros(1,length(po.(kpv{i}).slp));
    po.(kpv{i}).slpp=zeros(1,length(po.(kpv{i}).slp));

    %slope check - negative slope is to the south positive slope is to the
    %north
    for j=1:length(po.(kpv{i}).slp)
        if po.(kpv{i}).slp(j)>-7 && po.(kpv{i}).slp(j)<7
            po.(kpv{i}).slpln(j)=po.(kpv{i}).slp(j);
        elseif po.(kpv{i}).slp(j)<-7
            po.(kpv{i}).slpn(j)=po.(kpv{i}).slp(j); %rows that exceed south slope
        else
            po.(kpv{i}).slpp(j)=po.(kpv{i}).slp(j); %rows that exceed north slope
        end
    end
    if const.Nevados==1
        pbt = table(pmd.(kpv{i}).sect,pmd.(kpv{i}).row,po.(kpv{i}).bpx,po.(kpv{i}).bpy,po.(kpv{i}).upslope,po.(kpv{i}).downslope,...
            'VariableNames',["Section #","Row #","X Location","Y Location","Upslope (N+)","Downslope (N+)"]);
        pbt=rmmissing(pbt);
        if const.writefiles == 1
            if const.simple_out==0
                mkdir([const.outpath '/nevados' ]);
                writetable(pbt,[const.outpath '/nevados/' char(const.sect_name) '_slope_table.csv'])
            else
                mkdir([const.outpath '/nevados' ]);
                writetable(pbt,[const.outpath '/nevados/' char(const.sect_name) '_slope_table.csv'])
            end
        end
    end
    po.(kpv{i}).slpln(po.(kpv{i}).slpln==0)=NaN;
    po.(kpv{i}).slpn(po.(kpv{i}).slpn==0)=NaN;
    po.(kpv{i}).slpp(po.(kpv{i}).slpp==0)=NaN;
    po.(kpv{i}).os_plt=rmmissing([po.(kpv{i}).type, po.(kpv{i}).bpx, po.(kpv{i}).bpy, po.(kpv{i}).bpz, po.(kpv{i}).tpz, po.(kpv{i}).lpmax]);
    po.(kpv{i}).os_pls=rmmissing([po.(kpv{i}).type, po.(kpv{i}).bpx, po.(kpv{i}).bpy, po.(kpv{i}).bpz, po.(kpv{i}).tpz, po.(kpv{i}).lpmin]);
    po.(kpv{i}).os_pl=[po.(kpv{i}).os_plt; po.(kpv{i}).os_pls];
    po.(kpv{i}).os_ns=rmmissing([po.(kpv{i}).type, po.(kpv{i}).bpx, po.(kpv{i}).bpy, po.(kpv{i}).bpz, po.(kpv{i}).tpz, po.(kpv{i}).slpn']);
    po.(kpv{i}).os_ps=rmmissing([po.(kpv{i}).type, po.(kpv{i}).bpx, po.(kpv{i}).bpy, po.(kpv{i}).bpz, po.(kpv{i}).tpz, po.(kpv{i}).slpp']);
    po.(kpv{i}).oos=[po.(kpv{i}).os_pl; po.(kpv{i}).os_ns; po.(kpv{i}).os_ps];
    po.(kpv{i}).allpb=rmmissing([po.(kpv{i}).bpy,po.(kpv{i}).bpx,po.(kpv{i}).bpz]);
end
