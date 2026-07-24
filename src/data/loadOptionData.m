function data = loadOptionData(callsDir, putsDir)
%LOADOPTIONDATA Load raw call and put option CSV files.
%
% Output:
%   data(j).maturityCode
%   data(j).dates
%
%   data(j).callStrikes
%   data(j).calls
%
%   data(j).putStrikes
%   data(j).puts
%
% Important:
%   We do not require calls and puts to have the same strike grid.
%   Common strikes will be selected later only when needed.

    callFiles = dir(fullfile(callsDir, '*.csv'));
    putFiles  = dir(fullfile(putsDir, '*.csv'));

    callNames = string({callFiles.name});
    putNames  = string({putFiles.name});

    commonNames = intersect(callNames, putNames, 'stable');

    if isempty(commonNames)
        error('No common CSV files found between call and put folders.');
    end

    data = [];

    for j = 1:numel(commonNames)

        fileName = commonNames(j);

        callPath = fullfile(callsDir, fileName);
        putPath  = fullfile(putsDir, fileName);

        % Read call and put surfaces using your existing function
        callSurface = readOptionCsv(callPath, "Call");
        putSurface  = readOptionCsv(putPath, "Put");

        % Dates must match, because each row corresponds to a market date
        if ~isequal(callSurface.dates, putSurface.dates)
            error('Date grids do not match for file %s.', fileName);
        end

        % The maturity code should be the same because the file name is the same
        [~, maturityCode, ~] = fileparts(fileName);

        item = struct();

        item.maturityCode = string(maturityCode);
        item.dates = callSurface.dates;

        % Store calls separately
        item.callStrikes = callSurface.strikes(:);
        item.calls = callSurface.prices;

        % Store puts separately
        item.putStrikes = putSurface.strikes(:);
        item.puts = putSurface.prices;

        if isempty(data)
            data = item;
        else
            data(end+1) = item; %#ok<AGROW>
        end

    end

end