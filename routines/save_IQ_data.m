function save_IQ_data(IQData)
% function save_IQ_data(IQData)
% 
% Author: Yufeng Deng (Duke Unviersity)
% LICENSE: MIT

global filedir filetime comment

display('Start saving IQ data');

saveChannelData = evalin('base','saveChannelData');

tic

IQdata = squeeze(IQData);
I = real(IQdata);
Q = imag(IQdata);

if ~saveChannelData
    comment = '';
    filetime = datestr(clock, 'yyyymmddHHMMSS');
    
    CHMAT = fullfile(filedir,[filetime comment '_parameters']);
    Resource = evalin('base','Resource');
    PData = evalin('base','PData');
    Trans = evalin('base','Trans');
    TW = evalin('base','TW');
    TX = evalin('base','TX');
    Receive = evalin('base','Receive');
    ne = evalin('base','ne');
    nrefs = evalin('base','nrefs');
    T = evalin('base','T');
    T_idx = evalin('base','T_idx');
    lastBmodeEvent = evalin('base','lastBmodeEvent');
    lastBmodeReceive = evalin('base','lastBmodeReceive');
    lastBmodeTransmit = evalin('base','lastBmodeTransmit');
    pushAngleDegree = evalin('base','pushAngleDegree');
    Resource = rmfield(Resource,'DisplayWindow');
    Vvalue = 0;Vpush = 0;
    save(CHMAT,'Resource','PData','Trans','TW',...
        'TX','Receive','ne','nrefs','T','T_idx',...
        'lastBmodeEvent','lastBmodeReceive','lastBmodeTransmit',...
        'pushAngleDegree','Vvalue','Vpush','-v7.3')
end

IBIN = fullfile(filedir,[filetime comment '_IQreal.bin']);
QBIN = fullfile(filedir,[filetime comment '_IQimag.bin']);
    
fid=fopen(IBIN,'wb');
fwrite(fid,I,'int32');
fclose(fid);

fid=fopen(QBIN,'wb');
fwrite(fid,Q,'int32');
fclose(fid);

disp(['IQ data saved. Elapsed time is ' num2str(toc) ' seconds']);
% close all
end