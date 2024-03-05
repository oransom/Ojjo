function[drow, dpile, const] = a_4_npspans_6(const)

Folder = string(const.keypile);
if isfolder(Folder) ~= 1
    Message = sprintf('Error: The following folder does not exist:\n%s', Folder);
    uiwait(warndlg(Message));
    return;
end
filePattern = fullfile(Folder, '*.xlsx');
pileFiles   = dir(filePattern);
for k = 1:length(pileFiles)
    baseFileName = pileFiles(k).name;
    fullFileName = fullfile(Folder, baseFileName);
    kpv{k}=erase(pileFiles(k).name, ".xlsx");
    try
        dat.(kpv{k}) = readtable(fullFileName);
        dat.(kpv{k}).block(1:numel(dat.(kpv{k}).tpx))=string(kpv{k}); %insert the block name into the table dat
    catch ME
        if (strcmp(ME.identifier,'MATLAB:spreadsheet:book:fileOpen'))
            msg = 'Close all dependent project files';
            causeException = MException('MATLAB:spreadsheet:book:fileOpen',msg);
            ME = addCause(ME,causeException);
        end
        rethrow(ME)
    end
end
if sum(contains(fieldnames(dat.(kpv{1})),'row_length'))==0
    %% load in all span types
    Folder = string(const.spans);
    if isfolder(Folder) ~= 1
        Message = sprintf('Error: The following folder does not exist:\n%s', Folder);
        uiwait(warndlg(Message));
        return;
    end
    filePattern = fullfile(Folder, '*.xlsx');
    spanFiles   = dir(filePattern);
    for k = 1:length(spanFiles)
        baseFileName = spanFiles(k).name;
        fullFileName = fullfile(Folder, baseFileName);
        sv{k}=erase(spanFiles(k).name, ".xlsx");
        try
            span.(sv{k}) = readtable(fullFileName);
        catch ME
            if (strcmp(ME.identifier,'MATLAB:spreadsheet:book:fileOpen'))
                msg = 'Close all dependent project files';
                causeException = MException('MATLAB:spreadsheet:book:fileOpen',msg);
                ME = addCause(ME,causeException);
            end
            rethrow(ME)
        end
    end
else
    for i=1:length(kpv)
        [foo,fooi]=unique(dat.(kpv{i}).span);
        for j=1:length(foo)
            span.(foo{j}).s=foo{j};
            if span.(foo{j}).s<140
                span.(foo{j}).kps=[0:dat.(kpv{i}).row_length(fooi(j))/6:dat.(kpv{i}).row_length(fooi(j))];
                span.(foo{j}).kps=span.(foo{j}).kps';
                Alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
                for k=1:length(span.(foo{j}).kps)
                    span.(foo{j}).tl{k,1}=Alphabet(k);
                    if k~=4
                        span.(foo{j}).tt{k,1}='int';
                    else
                        span.(foo{j}).tt{k,1}='mtr';
                    end
                end
            elseif span.(foo{j}).s>140 && span.(foo{j}).s<250
                span.(foo{j}).kps=[0:dat.(kpv{i}).row_length(fooi(j))/10:dat.(kpv{i}).row_length(fooi(j))];
                span.(foo{j}).kps=span.(foo{j}).kps';
                Alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
                for k=1:length(span.(foo{j}).kps)
                    span.(foo{j}).tl{k,1}=Alphabet(k);
                    if k~=6
                        span.(foo{j}).tt{k,1}='int';
                    else
                        span.(foo{j}).tt{k,1}='mtr';
                    end
                end
            elseif span.(foo{j}).s>250 && span.(foo{j}).s<350
                span.(foo{j}).kps=[0:dat.(kpv{i}).row_length(fooi(j))/12:dat.(kpv{i}).row_length(fooi(j))];
                span.(foo{j}).kps=span.(foo{j}).kps';
                Alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
                for k=1:length(span.(foo{j}).kps)
                    span.(foo{j}).tl{k,1}=Alphabet(k);
                    if k~=7
                        span.(foo{j}).tt{k,1}='int';
                    else
                        span.(foo{j}).tt{k,1}='mtr';
                    end
                end
            else
                span.(foo{j}).kps=[0:dat.(kpv{i}).row_length(fooi(j))/14:dat.(kpv{i}).row_length(fooi(j))];
                span.(foo{j}).kps=span.(foo{j}).kps';
                Alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
                for k=1:length(span.(foo{j}).kps)
                    span.(foo{j}).tl{k,1}=Alphabet(k);
                    if k~=8
                        span.(foo{j}).tt{k,1}='int';
                    else
                        span.(foo{j}).tt{k,1}='mtr';
                    end
                end
            end
        end
    end
end

c = cellfun(@(fn) dat.(fn), fieldnames(dat), 'UniformOutput', false); %vertcat all NP files
drow = vertcat(c{:}); %dat_all is all NPS in vertical table
drow = renamevars(drow,["tpx","tpy"],["npx","npy"]);
for i=1:height(drow)
    drow.motorloc(i)=find(contains(span.(string(drow.span(i))).tt,'m','IgnoreCase',true));
    drow.mpy(i) = drow.npy(i)-span.(string(drow.span(i))).kps(drow.motorloc(i)); %any pile description with a m in it will be treated as a motor pile
    drow.spy(i) = drow.npy(i)-span.(string(drow.span(i))).kps(end); %southern pile y
    drow.mdptx(i)=drow.npx(i); %midpoint x
    drow.mdpty(i)=drow.npy(i)-span.(string(drow.span(i))).kps(end)/2; %midpoint y
end
const.xspacingpadded=const.ewdist+5; %average row spacing in x padded by 10 feet
if const.num_neigh==0
    const.xspacingpadded=const.xspacingpadded;
else
    const.xspacingpadded=(const.num_neigh/2)*(const.ewdist)+10;
end
rsx=[drow.mdptx,drow.mdpty]; %create array of x/y's to search against each other
rsy=[drow.mdptx,drow.mdpty];
[Idx,D]=rangesearch(rsx,rsy,const.xspacingpadded); %search for all row neighbors within padded x neighbor_distance
maxLengthCell=max(cellfun('size',Idx,2));  %finding the longest vector in the cell array
for m=1:length(Idx)
    for n=cellfun('size',Idx(m),2)+1:maxLengthCell
        Idx{m}(n)=0;   %zeropad the elements in each cell array with a length shorter than the maxlength
    end
end
neighbors=cell2mat(Idx); %turning IDX into numerical values
neighbors(neighbors==0)=NaN; %turning 0s into NaN's so we can avoid working with them
[m,n]=size(neighbors);
drowh.rnbrw=zeros(m,1); %preallocating arrays
drowh.rnbrc=1:1:m;
drowh.rnbrc=drowh.rnbrc';
drowh.rnbre=zeros(m,1);
if const.num_neigh==0
    for k=1:m
        for l=2:n
            if ~isnan(neighbors(k,l))
                if drow.mdptx(neighbors(k,l))<drow.mdptx(neighbors(k,1))...
                        && abs(drow.mdpty(neighbors(k,l))-drow.mdpty(neighbors(k,1)))<20 %finding neigbors to west
                    drowh.rnbrw(k,1)=neighbors(k,l);
                elseif drow.mdptx(neighbors(k,l))>drow.mdptx(neighbors(k,1))...
                        && abs(drow.mdpty(neighbors(k,l))-drow.mdpty(neighbors(k,1)))<20 % finding neighbors to east
                    drowh.rnbre(k,1)=neighbors(k,l);
                end
            end
        end
    end
else
    for k=1:m
        for l=2:n
            if ~isnan(neighbors(k,l))
                if drow.mdptx(neighbors(k,l))<drow.mdptx(neighbors(k,1))...
                        && abs(drow.mdpty(neighbors(k,l))-drow.mdpty(neighbors(k,1)))<20 %finding neigbors to west
                    drowh.rnbrw(k,l-1)=neighbors(k,l);
                elseif drow.mdptx(neighbors(k,l))>drow.mdptx(neighbors(k,1))...
                        && abs(drow.mdpty(neighbors(k,l))-drow.mdpty(neighbors(k,1)))<20 % finding neighbors to east
                    drowh.rnbre(k,l-1)=neighbors(k,l);
                end
            end
        end
    end
end
drowh.rnbrw(drowh.rnbrw==0)=NaN;
drowh.rnbre(drowh.rnbre==0)=NaN;
drow.row=drowh.rnbrc;
drow.nw_id=drowh.rnbrw;
drow.ne_id=drowh.rnbre;
drow.nnw=max(drowh.rnbrw,[],2);
drow.nne=min(drowh.rnbre,[],2);
%figure spacing out right here
nane=find(~isnan(drow.nne));
nanw=find(~isnan(drow.nnw));
spe=abs(drow.npx(nane)-drow.npx(drow.nne(nane)));
spw=abs(drow.npx(nanw)-drow.npx(drow.nnw(nanw)));
spacing=[spe;spw];
const.xspacing=mode(nonzeros(floor(spacing)));

%dpile = table["pnum","tpx","tpy","sect","row","pspan","block"];
si = 1;
for i=1:height(drow)
    ei=numel(span.(string(drow.span(i))).kps);
    ei=si+ei-1;
    dpile.bpx(si:ei,1)=drow.npx(i);
    dpile.bpy(si:ei,1)=drow.npy(i)-span.(string(drow.span(i))).kps;
    dpile.bpz(si:ei,1)=0;
    dpile.tpx(si:ei,1)=drow.npx(i);
    dpile.tpy(si:ei,1)=drow.npy(i)-span.(string(drow.span(i))).kps;
    dpile.tpz(si:ei,1)=0;
    dpile.sect(si:ei,1)=drow.sect(i);
    dpile.row(si:ei,1)=drow.row(i);
    dpile.pt(si:ei,1)=span.(string(drow.span(i))).tt;
    dpile.span(si:ei,1)=drow.span(i);
    dpile.block(si:ei,1)=drow.block(i);
    drow.si(i)=si;
    drow.ml(i)=drow.si(i)-1+drow.motorloc(i);
    drow.ei(i)=ei;
    si=ei+1;
end
dpile.bsr=strcat(dpile.block,'_',string(dpile.sect),'_',string(dpile.row)); %create unique identifier
dpile.pt=strrep(dpile.pt,'std','int'); %fixing other peoples mistakes
dpile=struct2table(dpile);