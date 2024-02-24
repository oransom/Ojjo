function sendOutlookMail(to,subject,body,attachments,cc)
% Sends email using MS Outlook.
%
% Modification of sendolmail to facilitate multiline emails.

% INPUTS:
% to
%    Recipient. Specified as a string or cell array of strings for multiple recipients.
%
% subject
%    A string indicating the subject of the email.
%
% body
%    The body of the email. If 'body' is a cell array of strings, this
%    version will parse the array and create a delimited string suitable
%    for display of multi-line (or multi-paragraph) emails. (This contrasts
%    with the original, which concatenated all input strings as a single
%    string.)
%
% attachments
%    A cell array of files which are to be attached. Use empty brackets if
%    you have no attachments, but you do have cc's. (E.g., {})
%
% cc
%    A string or cell array of strings indicating who should be cc'd.
%
% Note that the authorship of the original sendolmail seems to have been
% lost over the years. Was it Jeremy Nersasian (who shared the function
% when he was in Technical Support at MathWorks)?
%
% Brett Shoelson, PhD
% brett.shoelson@mathworks.com
% 12/04/2018

if nargin < 3
    error('Function %s requires at least three input arguments.\n', mfilename);
end
% Create object and set parameters.
h = actxserver('outlook.Application');
mail = h.CreateItem('olMail');
mail.Subject = subject;
if isa(to,'cell')
    to = strjoin(to,';');
end
mail.To = to;
mail.BodyFormat = 'olFormatHTML';
%'<p>First line here</p><p>Second line here</p>'
if isa(body,'cell')
    body = reformatBody(body);
end
mail.HTMLBody = body;

% Add attachments, if specified.
if nargin > 3 && ~isempty(attachments)
    for ii = 1:length(attachments)
        mail.attachments.Add(attachments{ii});
    end
end
if nargin > 4 
    if isa(cc,'cell')
        cc = strjoin(cc,';');
    end
else 
    cc = [];
end
if ~isempty(cc)
    mail.cc = cc;
end

% Send message and release object.
mail.Send;
h.release;
end

function body = reformatBody(body)
body = sprintf('<p>%s</p>',body{:});
end %reformatBody