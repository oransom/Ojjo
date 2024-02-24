clear all
close all
clc
setpref('Internet','E_mail','natsilane@orcas.sh');
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username','natsilane@orcas.sh');
setpref('Internet','SMTP_Password','ylvbaodsnvanxkbr');
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');


%nevados
%directoryPath = "D:\ORCaS - Drive\My Drive\Nevados\Auto Run In"; % Replace with your directory path
%runPath = "D:\ORCaS - Drive\My Drive\Nevados\Auto Run Out";
%RPCS
directoryPath = "D:\RPCS_Box\Box Sync\Planar\1A) Auto_Input_Windows"; % Replace with your directory path
runPath = "D:\RPCS_Box\Box Sync\Planar\1B) Auto_Output";
allFilesOld = dir(fullfile(directoryPath, '*.zip')); % This will get a list of all zip files in the directory
if ~isempty(allFilesOld)
    for i=1:length(allFilesOld)
        aFO_namelist{i}=allFilesOld(i).name;
    end
else
    aFO_namelist=[];
end
hasBeenRun = allFilesOld;
aFN_namelist={};
lengthLastMsg=fprintf('Starting monitoring \n');
while true
    try
        ctime=timeofday(datetime('now'));
        starttime='08:00:00'; endtime='18:00:00';
        if ctime > '09:00:00' & ctime <'09:05:01'...
            | ctime > '10:00:00' & ctime <'10:05:01'...
            | ctime > '11:00:00' & ctime <'11:05:01'...
            | ctime > '12:00:00' & ctime <'12:05:01'...
            | ctime > '13:00:00' & ctime <'13:05:01'...
            | ctime > '14:00:00' & ctime <'14:05:01'...
            | ctime > '15:00:00' & ctime <'15:05:01'...
            | ctime > '16:00:00' & ctime <'16:05:01'...
            | ctime > '17:00:00' & ctime <'17:05:01'...
            | ctime > '18:00:00' & ctime <'18:05:01'...
            | ctime > '19:00:00' & ctime <'19:05:01'...
            | ctime > '20:00:00' & ctime <'20:05:01'...
            | ctime > '21:00:00' & ctime <'21:05:01'...
            | ctime > '23:00:00' & ctime <'23:59:59'
            %movefolders(directoryPath);
            pause(600)
        elseif ctime > starttime & ctime < endtime
            pause(600)
        else
            pause(3600)
        end
        allFilesNew = dir(fullfile(directoryPath, '*.zip')); % This will get a list of all zip files in the directory
        if ~isequal(allFilesOld,allFilesNew)
            lengthLastMsg=fprintf('new files found\n');
            pause(.5)
            ii=0;
            for i = 1:length(allFilesNew)
                aFN_namelist{i}=allFilesNew(i).name;
                if ~ismember(aFN_namelist{i},aFO_namelist)
                    fprintf('Unzipping and running folder %d: %s\n', i, allFilesNew(i).name);
                    ii=ii+1;
                    newrun{ii}=horzcat(allFilesNew(i).folder,'/',allFilesNew(i).name);
                    for j=1:length(newrun)
                        unzip(newrun{j},runPath)
                        pause(1)
                        Folder=string(strcat(runPath,'/',extractBefore(allFilesNew(i).name,'.zip')));
                        filePattern = fullfile(Folder, '*.xlsx');
                        cf   = dir(filePattern);
                        for k=1:length(cf)
                            if ~startsWith(cf(k).name,'~$')
                                config_file=horzcat(cf(k).folder,'/',cf(k).name);
                                aauto_MainProgram(config_file);
                                myaddress='owen@orcas.sh';
                                subject='RPCS Monitoring Program Finished';
                                message=['RPCS Planar Monitoring Program has run on: ', allFilesNew(i).name];
                                sendmail(myaddress, subject, message);
                            end
                        end
                    end
                end
            end
            aFO_namelist=aFN_namelist;
            allFilesOld=allFilesNew;
            clear newrun
            lengthLastMsg=fprintf('\n all new runs completed at %s\n',timeofday(datetime('now')));
        else
            fprintf(repmat('\b', 1, lengthLastMsg));
            lengthLastMsg=fprintf('no new files %s\n',timeofday(datetime('now')));
        end
    catch ME
        fprintf('Error in monitoring program... resetting');
        myaddress='owen@orcas.sh';
        subject='RPCS Monitoring Program Error';
        message=['RPCS Planar Monitoring Program has errored out error given as: ', ME.message];
        %sendmail(myaddress, subject, message);
        continue
    end
end

function [] = movefolders(directoryPath)
dir = directoryPath;
tdir='D:\RPCS_Box\Box Sync\Planar\Auto_Input_Archive';
cdir=pwd;
cd(dir)
foo=datestr(today);
movefile ([dir,'\*'], [tdir, '\', foo, '\']);
cd(cdir)
end
