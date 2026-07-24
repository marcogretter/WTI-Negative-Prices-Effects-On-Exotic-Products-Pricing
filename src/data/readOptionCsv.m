function surface = readOptionCsv(filePath, optionType)
%READOPTIONCSV Read one option CSV file.
% Expected CSV format:
%   first row    = strike headers
%   first column = dates
%   body         = option prices
%
% Output struct:
%   maturityCode, optionType, dates, strikes, prices

    raw = readcell(filePath);

    if size(raw,1) < 2 || size(raw,2) < 2
        error('Invalid CSV structure: %s', filePath);
    end

    [~, fileName, ~] = fileparts(filePath);

    strikes = cellfun(@toDouble, raw(1, 2:end));
    dates   = parseDates(raw(2:end, 1));
    prices  = cellfun(@toDouble, raw(2:end, 2:end));

    validStrikes = isfinite(strikes);

    surface = struct();
    surface.maturityCode = string(fileName);
    surface.optionType   = string(optionType);
    surface.filePath     = string(filePath);
    surface.dates        = dates(:);
    surface.strikes      = strikes(validStrikes);
    surface.prices       = prices(:, validStrikes);

end

function x = toDouble(v)
    if isempty(v)
        x = NaN;
    elseif isnumeric(v)
        x = double(v);
    else
        x = str2double(strtrim(string(v)));
    end
end

function dates = parseDates(dateCells)
    n = numel(dateCells);
    dates = NaT(n,1);

    for i = 1:n
        v = dateCells{i};

        if isdatetime(v)
            dates(i) = v;

        elseif isnumeric(v)
            % The dataset stores dates as yyyymmdd, e.g. 20200106.
            if v > 1e7
                dates(i) = datetime(string(sprintf('%.0f', v)), ...
                    'InputFormat', 'yyyyMMdd');
            else
                dates(i) = datetime(v, 'ConvertFrom', 'excel');
            end

        else
            s = strtrim(string(v));

            if strlength(s) == 8 && all(isstrprop(char(s), 'digit'))
                dates(i) = datetime(s, 'InputFormat', 'yyyyMMdd');
            else
                try
                    dates(i) = datetime(s, 'InputFormat', 'yyyy-MM-dd');
                catch
                    dates(i) = datetime(s);
                end
            end
        end
    end

    if any(isnat(dates))
        error('Some dates could not be parsed.');
    end
end