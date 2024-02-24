function emailer2(const, kpv, pdfout, mname, pd_out, rp_loc, pbt)

for i=1:length(kpv)
    sslope=unique(pd_out.(kpv{i}).slp_deg(pd_out.(kpv{i}).slp_deg<0));
    nslope=unique(pd_out.(kpv{i}).slp_deg(pd_out.(kpv{i}).slp_deg>0));
end

sspct=length(sslope)/(length(sslope)+length(nslope));
nspct=length(nslope)/(length(sslope)+length(nslope));

sendto=string(const.email{1}.sendto(1));
hcc=string(const.email{1}.sendto(1));
if length(const.email{1}.sendto)>1
    cc=const.email{1}.sendto(2:end)';
    hcc=const.email{1}.sendto(1:end)';
    d = [hcc',[repmat({','},numel(hcc)-1,1);{[]}]]';
    e = [d{:}];
end

subjecth=['Planar and associated information for: ', const.project, ' ', const.tracker];

pbtc=table2cell(pbt);
pbth=pbt.Properties.VariableNames;
Matrix1 = vertcat(pbth,pbtc);%
% rand(10,5,2); % 3D numerical
out = HTMLtable(Matrix1(:,:));%,'ShowOutput',true);
ems = HTMLtable((cellstr(hcc))');
formatspec="<p>This email has the following receipients: <p><p>%s<p><p> Here is the %s %s planar you requested. The max north-south slope" + ...
    " for the project is %4.2f degrees with a %4.2f to %4.2f ratio of south to north facing trackers. <p>" + ...
    " The percentage of each TTH is shown in the table below.<p><p> %s <p> Please let us know if you have any questions <p><p>" + ...
    "www.orcas.sh (__-){ <p>" + ...
    "+1.415.450.7558 <p>";
message1=sprintf(formatspec,string(ems),string(const.project),string(const.tracker),max(abs(pd_out.(kpv{i}).slp_deg)),sspct*100,nspct*100,string(out));


if const.Nevados ~=1 && ~contains(const.tracker,'ojjo','IgnoreCase',true)
    attachments={pdfout,mname,rp_loc};
elseif contains(const.tracker,'ojjo','IgnoreCase',true)
    attachments={pdfout};
else
    attachments={pdfout,mname};
end
subjecth2=[subjecth{:}];
subject=string(subjecth2);
%% this is the outlook version
%sendmail(sendto, subject, message, attachments);
% if length(const.email{1}.sendto)>1
%     sendOutlookMail(sendto,subject,message,attachments,cc)
% else
%     sendOutlookMail(sendto,subject,message,attachments)
% end
%% this is the new os agnostic html version

setpref('Internet','E_mail','natsilane@orcas.sh');
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username','natsilane@orcas.sh');
setpref('Internet','SMTP_Password','ylvbaodsnvanxkbr');
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');

% server='smtp.gmail.com';
% setpref('Internet','E_mail','natsilane@orcas.sh');
% setpref('Internet','SMTP_Server','smtp.gmail.com');
% setpref('Internet','SMTP_Username','natsilane@orcas.sh');
% setpref('Internet','SMTP_Password','ylvbaodsnvanxkbr');
% props = java.lang.System.getProperties;
% props.setProperty( 'mail.smtp.auth', 'true' );
% props.setProperty( 'mail.smtp.user', 'natsilane@orcas.sh' );
% props.setProperty( 'mail.smtp.password', 'ylvbaodsnvanxkbr');
% props.setProperty( 'mail.smtp.host', 'smtp.gmail.com');
% props.setProperty( 'mail.smtp.port', '587' );
% props.setProperty( 'mail.smtp.starttls.enable', 'true' );

for i=1:length(hcc)
    tofoo=hcc{i};
    sendhtmlemail(tofoo,subject,message1,attachments);
end


