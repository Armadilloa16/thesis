function sortedPeaks = sortPeaks(L,LXY) 

    if nargin < 3
        includeVars = [];
    end
    % Compile all peaks into one list
    allPeaks = zeros(sum(LXY(:,3)),3+length(includeVars));
    curLoc = 1;
    for spec = 1:size(LXY,1)
        allPeaks(curLoc:(curLoc+LXY(spec,3)-1),1) = spec*ones(LXY(spec,3),1);
        allPeaks(curLoc:(curLoc+LXY(spec,3)-1),2) = 1:LXY(spec,3);
        allPeaks(curLoc:(curLoc+LXY(spec,3)-1),3) = L{spec}(:,1);
        for i = 1:length(includeVars)
            allPeaks(curLoc:(curLoc+LXY(spec,3)-1),3+i) = L{spec}(:,includeVars(i));
        end
        curLoc = curLoc + LXY(spec,3);
    end

    % Sort allPeaks by m/z location into sortedPeaks
    [~,I] = sort(allPeaks(:,3));
    sortedPeaks = allPeaks(I,:);
    
end