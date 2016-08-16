function [mcdata,vBinCentrs,mdataL] = bin1dFast(sortedPeaks,binSize,offCAdj)
    
    sortedSpectra = sortedPeaks(:,1);
    genMdataL = nargout > 2;
    if genMdataL
        sortedPeakNo = sortedPeaks(:,2); 
    end
    sortedPeaks = sortedPeaks(:,3);

    %% Check inputs
    % Check nargin, binSize, and offCAdj
    if nargin < 2
        error('ERROR: !!! Not enough arguments for bin1dVec')
    elseif nargin >= 2
        if ~isnumeric(binSize)
            error('ERROR: !!! Bin size not numeric')
        elseif size(binSize,1) + size(binSize,2) ~= 2
            error('ERROR: !!! tried to select multiple bin sizes - not currently supported')
        elseif nargin == 2
            offCAdj = 0;
        elseif nargin == 3
            if ~isnumeric(offCAdj)
                error('ERROR: !!! Wiggle parameter not numeric')
            elseif size(offCAdj,1) + size(offCAdj,2) ~= 2
                error('ERROR: !!! tried to select multiple wiggle parameters - not currently supported')
            elseif offCAdj < -1*(binSize/2) 
                error('ERROR: !!! offCAdj argument to bin1dVec is negative - not allowed')
            elseif offCAdj > (binSize/2)
                error('ERROR: !!! offCAdj argument  is larger than or equal to binSize argument to bin1dVec - not allowed')
            end
        elseif nargin > 3
            error('ERROR: !!! Too many arguments for bin1dVec')
        end
    end
    
    % Check that sortedPeaks is sorted
%     [sortedPeaks,I] = sort(sortedPeaks);
%     sortedSpectra = sortedSpectra(I);
    
    %% Do stuff
    % Wiggle peaks
    sortedPeaks = sortedPeaks - offCAdj;
    % Find range
    minBinCentr = correctRounding(floor((sortedPeaks(1))/binSize)*binSize);
    maxBinCentr = correctRounding(ceil(sortedPeaks(length(sortedPeaks))/binSize)*binSize);
    vBinCentrs = correctRounding(minBinCentr:binSize:maxBinCentr);
    
    % Create vector of binary data - non-emptry m/z bins
    vbdata = false(length(vBinCentrs),1);
    for peakIdx = 1:length(sortedPeaks)
        binC = correctRounding(round(sortedPeaks(peakIdx)/binSize - 10^-10)*binSize);
        binIdx = round(((binC - minBinCentr)/binSize) + 1);
        vbdata(binIdx) = true(1,1);
    end
    
    % Create matrix of count data - non-empty m/z bins
    vBinCentrs = vBinCentrs(vbdata);
    mcdata = zeros(sum(vbdata),length(unique(sortedSpectra)));
    if genMdataL 
    mdataL = mcdata;
    end
    for peakIdx = 1:length(sortedPeaks)
        binC = correctRounding(round(sortedPeaks(peakIdx)/binSize - 10^-10)*binSize);
        binIdx = find(vBinCentrs == binC);
        mcdata(binIdx,sortedSpectra(peakIdx)) = mcdata(binIdx,sortedSpectra(peakIdx)) + 1;
        if genMdataL
            if mdataL(binIdx,sortedSpectra(peakIdx)) == 0
                mdataL(binIdx,sortedSpectra(peakIdx)) = sortedPeakNo(peakIdx);
            end
        end
    end
    
    vBinCentrs = correctRounding(vBinCentrs + offCAdj);
    
end