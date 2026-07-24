function mktData = extract_data(projectRoot,valuationDate)
% Load raw option data and futures expiry calendar for the valuation date.
%
% Returns mktData with fields:
%   snapshot        struct array, one element per maturity
%   futureExpiries  datetime column vector of all maturities
%   futureCodes     string column vector of maturity codes
%   valuationDate   datetime scalar

callsDir = fullfile(projectRoot, 'Data', 'datacalls');
putsDir  = fullfile(projectRoot, 'Data', 'dataputs');

%valuationDate = datetime(2020, 6, 2);

% Load raw option data
data = loadOptionData(callsDir, putsDir);

% Extract valuation-date snapshot and remove only NaN/Inf quotes
snapshot = getOptionSnapshot(data, valuationDate);

% Load futures expiries
expiriesFile = fullfile(projectRoot, 'Data', 'Expiries_Futures.txt');

fid = fopen(expiriesFile, 'r');
tmp = textscan(fid, '%s %s');
fclose(fid);

futureCodes = string(tmp{1});
futureExpiries = datetime(string(tmp{2}), 'InputFormat', 'yyyy/MM/dd');

mktData.snapshot       = snapshot;
mktData.futureExpiries = futureExpiries(:);
mktData.futureCodes    = futureCodes(:);
mktData.valuationDate  = valuationDate;
end