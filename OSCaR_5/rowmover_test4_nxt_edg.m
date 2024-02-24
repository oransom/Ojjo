function [rnbr_table, mover_table, fliped] = rowmover_test4_nxt_edg(const, ot)
% constant declarations
% const.row_delta=1.375;
neighbor_distance=const.neighbor_distance;
ns_off=20;
ewpad=0.1;
int2edg=0.065;
int2ext=0.12;
%int to edge - 6.5% - 3.72 deg
%int to ext - 12.0% - 6.84 deg

rsx=[ot.row_tpx_avg, ot.row_tpy_avg];
rsy=[ot.row_tpx_avg, ot.row_tpy_avg];
[Idx,D]=rangesearch(rsx,rsy,neighbor_distance); %search for all row neighbors within neighbor_distance
maxLengthCell=max(cellfun('size',Idx,2));  %finding the longest vector in the cell array
for i=1:length(Idx)
   for j=cellfun('size',Idx(i),2)+1:maxLengthCell
       Idx{i}(j)=0;   %zeropad the elements in each cell array with a length shorter than the maxlength
   end
end
neighbors=cell2mat(Idx); %turning IDX into numerical values
neighbors(neighbors==0)=NaN; %turning 0s into NaN's so we can avoid working with them

[m,n]=size(neighbors);
rnbr.tpxc=ot.row_tpx_avg;   rnbr.tpyc=ot.row_tpy_avg; rnbr.tpzc=ot.row_tpz_avg; rnbr.tpzo=ot.row_tpz_avg;
rnbr.wpAvg=ot.preveal_mean; rnbr.wpRemain=ot.preveal_remain;     rnbr.row_length=ot.row_length;
rnbr.w=zeros(m,1);          rnbr.e=zeros(m,1);        rnbr.bsrw=cell(m,1);     rnbr.bsrc=ot.bsr;
rnbr.bsre=cell(m,1);        rnbr.movecr=zeros(m,1);   rnbr.canmove=zeros(m,1);       %allocate matrices for row neighbor information

for i=1:m
    for j=2:n
        if ~isnan(neighbors(i,j))
            if rnbr.tpxc(neighbors(i,j))<rnbr.tpxc(neighbors(i,1)) && abs(rnbr.tpyc(neighbors(i,j))-rnbr.tpyc(neighbors(i,1)))<ns_off %finding neigbors to west
                rnbr.w(i)=neighbors(i,j);
                rnbr.bsrw(i)=ot.bsr(rnbr.w(i));
            elseif rnbr.tpxc(neighbors(i,j))>rnbr.tpxc(neighbors(i,1)) && abs(rnbr.tpyc(neighbors(i,j))-rnbr.tpyc(neighbors(i,1)))<ns_off % finding neighbors to east
                rnbr.e(i)=neighbors(i,j);
                rnbr.bsre(i)=ot.bsr(rnbr.e(i));
            end
        end
    end
end

rnbr.e(rnbr.e==0)=NaN; rnbr.w(rnbr.w==0)=NaN;   %Turn zeros to NaN's in neighbors
rnbr.d_west=zeros(m,1); rnbr.d_east=zeros(m,1);
rnbr.m_west=zeros(m,1); rnbr.m_east=zeros(m,1);
rnbr.slp_e=zeros(m,1); rnbr.slp_ep=zeros(m,1);
rnbr.slp_w=zeros(m,1); rnbr.slp_wp=zeros(m,1);
rnbr.flipexterior=zeros(length(rnbr.e),1);
rnbr.dw_ed_al=zeros(m,1); rnbr.dw_ex_al=zeros(m,1);
rnbr.de_ed_al=zeros(m,1); rnbr.de_ex_al=zeros(m,1);

for i=1:length(rnbr.tpzc)
    if ~isnan(rnbr.e(i))
        rnbr.slp_e(i)=atand((rnbr.tpzc(i)-rnbr.tpzc(rnbr.e(i)))/abs(rnbr.tpxc(i)-rnbr.tpxc(rnbr.e(i))));
        rnbr.slp_ep(i)=100*(rnbr.tpzc(i)-rnbr.tpzc(rnbr.e(i)))/abs(rnbr.tpxc(i)-rnbr.tpxc(rnbr.e(i)));
        rnbr.de_ed_al(i)=abs(rnbr.tpxc(i)-rnbr.tpxc(rnbr.e(i)))*int2edg; %delta edge allowable
        rnbr.de_ex_al(i)=abs(rnbr.tpxc(i)-rnbr.tpxc(rnbr.e(i)))*int2ext; %delta ext allowable
    end
    if ~isnan(rnbr.w(i))
        rnbr.slp_w(i)=atand((rnbr.tpzc(rnbr.w(i))-rnbr.tpzc(i))/abs(rnbr.tpxc(rnbr.w(i))-rnbr.tpxc(i)));
        rnbr.slp_wp(i)=100*(rnbr.tpzc(rnbr.w(i))-rnbr.tpzc(i))/abs(rnbr.tpxc(rnbr.w(i))-rnbr.tpxc(i));
        rnbr.dw_ed_al(i)=abs(rnbr.tpxc(rnbr.w(i))-rnbr.tpxc(i))*int2edg; %delta edge allowable
        rnbr.dw_ex_al(i)=abs(rnbr.tpxc(rnbr.w(i))-rnbr.tpxc(i))*int2ext; %delta ext allowable
    end
end

morecanmove=1; %So the loop will run once
rnbr.wpRemain_n = rnbr.wpRemain;
while morecanmove>0
    rnbr.canmove=rnbr.canmove*0;
    rnbr.movecr=rnbr.movecr*0;
for i=1:length(rnbr.w)
    if ~isnan(rnbr.w(i))
    if rnbr.wpRemain(i)>0
     if rnbr.tpzc(i)-rnbr.tpzc(rnbr.w(i))<-rnbr.dw_ed_al(i) %if tpz of crow is greater than const.row_delta of west neighbor row
         rnbr.movecr(i)=1; %binary 1/0 for out of spec
         rnbr.d_west(i)=rnbr.tpzc(i)-rnbr.tpzc(rnbr.w(i)); %delta between rows that violate
         rnbr.m_west(i)=abs(rnbr.d_west(i))-rnbr.dw_ed_al(i);
         if rnbr.m_west(i)<rnbr.wpRemain_n(i) %check to see if delta is greater than wp, if not, can move
            rnbr.canmove(i)=1;
         end
     end
    end
    end
end

for i=1:length(rnbr.e)
    if ~isnan(rnbr.e(i))
    if rnbr.wpRemain(i)>0
     if rnbr.tpzc(i)-rnbr.tpzc(rnbr.e(i))<-rnbr.de_ed_al(i) %if tpz of crow is greater than const.row_delta of east neighbor row
         rnbr.movecr(i)=1;
         rnbr.d_east(i)=rnbr.tpzc(i)-rnbr.tpzc(rnbr.e(i));
         rnbr.m_east(i)=abs(rnbr.d_east(i))-rnbr.de_ed_al(i);
         if rnbr.m_east(i)<rnbr.wpRemain_n(i) %check to see if delta is greater than wp, if not, can move
            rnbr.canmove(i)=1;
         end
     end
    end
    end
end

for i=1:length(rnbr.tpzc)
    if rnbr.canmove(i)==1
        if rnbr.m_west(i)>0 && rnbr.m_east(i)>0
            if rnbr.m_west(i)>rnbr.m_east(i)
                rnbr.tpzc(i)=rnbr.tpzc(i)+rnbr.m_west(i);
            else
                rnbr.tpzc(i)=rnbr.tpzc(i)+rnbr.m_east(i);
            end
        elseif rnbr.m_west(i)>0
            rnbr.tpzc(i)=rnbr.tpzc(i)+rnbr.m_west(i);
        elseif rnbr.m_east(i)>0
            rnbr.tpzc(i)=rnbr.tpzc(i)+rnbr.m_east(i);
        else
            rnbr.tpzc(i)=rnbr.tpzc(i);
        end
    end
end

for i=1:length(rnbr.tpzc)
    if rnbr.canmove(i)==1
        rnbr.wpRemain_n(i)=rnbr.wpRemain_n(i)-(rnbr.tpzc(i)-rnbr.tpzo(i));
    end
end
morecanmove=sum(rnbr.canmove);
movecrsum=sum(rnbr.movecr);
end

for i=1:length(rnbr.movecr)
    if rnbr.movecr(i)==1 && rnbr.d_west(i) < -rnbr.dw_ed_al(i) %check if current row is triggering and then flip the next row over
        if ~isnan(rnbr.w(i))
        rnbr.flipexterior(rnbr.w(i))=1;
        end
        if ~isnan(rnbr.w(rnbr.w(i)))
        if rnbr.tpzc(rnbr.w(rnbr.w(i)))-rnbr.tpzc(i)>rnbr.dw_ed_al(i) %look two rows over to see if we need to flip it as well based on tpzc there vs. center
            rnbr.flipexterior(rnbr.w(rnbr.w(i)))=1;
        end
        end
    elseif rnbr.movecr(i)==1 && rnbr.d_east(i) < -rnbr.de_ed_al(i)
        if ~isnan(rnbr.e(i))
        rnbr.flipexterior(rnbr.e(i))=1;
        end
        if ~isnan(rnbr.e(rnbr.e(i)))
        if rnbr.tpzc(rnbr.e(rnbr.e(i)))-rnbr.tpzc(i)>rnbr.de_ed_al(i)
            rnbr.flipexterior(rnbr.e(rnbr.e(i)))=1;
        end
        end
    end
end

fliped.bsr=rnbr.bsrc; fliped.tpxc=rnbr.tpxc; fliped.tpyc=rnbr.tpyc; fliped.tpzc=rnbr.tpzc;
fliped.row_len=rnbr.row_length;
fliped.flip=rnbr.flipexterior;
fliped.flip(fliped.flip==0)=NaN;
fliped=struct2table(fliped);
fliped=rmmissing(fliped);

% %ew slope histogram
% slopedist=rnbr.slp_e;
% figure ('Visible', 'off')
% set(gcf, 'Position',  [100, 100, 600, 400])
% edges=[floor(min(slopedist)):.25:ceil(max(slopedist))];
% edges2=[floor(min(slopedist))-1:1:ceil(max(slopedist))+1];
% [pbin pbedges]=histcounts(slopedist,edges2);
% pbinp=100*(pbin/(sum(pbin(:))));
% formatSpec = '%.1f';
% for j=1:size(pbinp,2)
%     pbins{j} = num2str(pbinp(j),formatSpec);
% end
% %pbt = table(e2name,pbins','VariableNames',["Pile Bin","Approximate Percentage"]);
% h1=histogram(slopedist,edges,'Normalization','probability'); %error here might mean you need to up your neighbor distance
% yl=(get(gca,'YLim'));
% set(gca,'YLim',[yl(1),yl(2)+.13*yl(2)]);
% ylabel('Percentage of Rows')
% xlabel('Row Slope - West Positive')
% n=get(h1,'Values');
% xn=get(h1,"BinEdges");
% formatSpec = '%.2f';
% barstrings=num2str(n',formatSpec);
% text(xn(2:end), n, barstrings,'horizontalalignment','center','verticalalignment','bottom','Rotation',70)
% %labels=num2str([[edges(1:end-1)].' [edges(2:end)].'],'%0.2f-%0.2f');
% % set labels
% %set(gca,'XTickLabel',labels,'XTickLabelRotation',30)
% 
% set(gca, 'units', 'normalized'); %Just making sure it's normalized
% Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
% %[Left Bottom Right Top] spacing
% NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
% set(gca, 'Position', NewPos);
% 
% if const.writefiles == 1
%     if const.simple_out==0
%         mkdir([const.outpath '/figures' ]);
%         saveas(gcf,[const.outpath '/figures/' char(const.sect_name) ' ew_hist.png'])
%     else
%         mkdir([const.outpath '/figures' ]);
%         saveas(gcf,[const.outpath '/figures/' char(const.sect_name) ' ew_hist.png'])
%     end
%     foo=1;
% end

if const.writefiles==1 && const.simple_out==1
    mkdir([const.outpath '/fliped' ]);
    writetable(fliped,[const.outpath '/fliped/' char(const.sect_name) ' fliped.csv'])
end

rnbr.raiserow=rnbr.tpzc-rnbr.tpzo;
mover.bsr=rnbr.bsrc;
mover.raiserow=rnbr.raiserow;
mover.raiserow(mover.raiserow==0)=NaN;
rnbr_table=struct2table(rnbr);
mover_table=struct2table(mover);
mover_table=rmmissing(mover_table);

if const.writefiles==1 && const.simple_out==1
    mkdir([const.outpath '/rowdata' ]);
    writetable(rnbr_table,[const.outpath '/rowdata/' char(const.sect_name) ' rowdata.csv'])
end

end
