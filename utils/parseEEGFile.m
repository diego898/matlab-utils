function [datesTimesEventCodes, eeg] = parseEEGFile(fileName, numChannels, samplingRate)
%parseEEGFile   Read in and parse EEG File
%
%   parseEEGFile(filename, numChannels) will read in the specified EEG File
%   which was exported from MJs EGG Machine.
%
%   Default values for numChannels = 33, samplingRate = 512
%
%   datesTimesEventCodes is a matrix of size numOfSamples * 7 where each
%   column is month-date-year-hour-minute-second-event_code
%
%   eeg is a matrix of size numOfSamples * numChannels
%
%   NOTE: SHORT Values APPEAR AS -INF, and the 1st second of recording is
%   thrown away
%
%   For formatting information view source



%% Notes/Settings
% The following line appears and is probably due to an error in the mach.
% ---,BREAK,IN,DATA,---,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

% SHORT appears in channels sometimes, inserting a -inf in its place

% Format of file by column:
% 1 - Date (String)
% 2 - Time (String)
% 3 - Event ID (Float)
% [4:37] - 33 Channels (Float)
% [38:?] - Trigger - Ignored till end of line

% The following options were used to parse the file:
% CollectOutput, true - prevents crazy cells from forming
% commentStyle, '---' - this is to ignore those BROKEN DATA strings 
% EmptyValue, -Inf - this replaces empty values with this so we can tell
% TreatAsEmpty, 'SHORT' - whenever you see 'SHORT', treat it as empty
% HeaderLines, 15 - Removed first 15 lines of header text



%% Default Values
switch nargin
    case 2
        samplingRate = 512;
    case 1
        numChannels = 33;
        samplingRate = 512;
    case 0
        error('You have to atleast provide a filename');
end



%% Parsing
fid = fopen(fileName,'r');
format = ['%u/%u/%u ' '%u:%u:%u ' '%u ' repmat('%f ', [1 numChannels]) '%*[^\n\r]'];
fileData = textscan(fid, format, 'CollectOutput', true, 'commentStyle', '---', 'EmptyValue', -Inf, 'TreatAsEmpty', 'SHORT', 'HeaderLines', 15);
fclose(fid);
datesTimesEventCodes = fileData{1};
eeg = fileData{2};

% Throw away first second since it might not contain all 512 samples
numRowsDisc = sum(datesTimesEventCodes(1:samplingRate,6) == datesTimesEventCodes(1,6));
datesTimesEventCodes(1:numRowsDisc, :) = [];
eeg(1:numRowsDisc, :) = [];



end     % end function