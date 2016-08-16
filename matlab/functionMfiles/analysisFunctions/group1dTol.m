function [mcdata,vBinCentrs,vBinWidths,mdataL] = group1dTol(sortedPeaks,tol)
    
    sortedSpectra = sortedPeaks(:,1);
    genMdataL = nargout > 3;
    if genMdataL
        sortedPeakNo = sortedPeaks(:,2); 
    end
    sortedPeaks = sortedPeaks(:,3);

    vdif = [find(diff(sortedPeaks)>tol) length(sortedPeaks)];
    if length(unique(sortedSpectra)) ~= max(sortedSpectra)
        error('Something went wrong with sortedSpectra')
    else
        nSpec = max(sortedSpectra);
    end
    mcdata = zeros(length(vdif),nSpec);
    vBinCentrs = zeros(length(vdif),1);
    vBinWidths = zeros(length(vdif),1);
    % This is slow, and could be recoded better, as bin1dFast has been for example.
    curLoc = 1;
    for i = 1:length(vdif)
        group = curLoc:vdif(i);
        curLoc = vdif(i)+1;
        minMZ = min(sortedSpectra(group));
        maxMZ = max(sortedSpectra(group));
        vBinCentrs(i) = mean([minMZ maxMZ]);
        vBinWidths(i) = maxMZ - minMZ;
        
        values = unique(
        for spec = 1:max(sortedSpectra)
            mcdata(i,spec) = sum(sortedPeaks(group + (sortedSpectra==spec)==2));
        end
    end
    
    
end