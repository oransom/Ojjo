
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
            movefolders;
            pause(30)
        elseif ctime > starttime & ctime < endtime
            pause(30)
        else
            pause(30)
        end

    end
end

function [] = movefolders(directoryPath)
        dir = 'E:\Dropbox (Personal)\ORCaS\ORCaS Code\Testbed\ORCaS Projects\testfolder2';
        tdir='E:\Dropbox (Personal)\ORCaS\ORCaS Code\Testbed\ORCaS Projects\testfolder1';
% dir = directoryPath;
% tdir='D:\RPCS_Box\Box Sync\Planar\Auto_Input_Archive';
cdir=pwd;
cd(dir)
foo=datestr(today);
movefile ([dir,'\*'], [tdir, '\', foo, '\']);
cd(cdir)
end