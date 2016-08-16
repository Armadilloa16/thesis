function vC = combineDims(vCin)
    %%% Takes input vC1, vC2, ... and combines them into vC
    
    %%% Check that the input is in the supported format.
    if size(vCin,1)+size(vCin,2)-(size(vCin,1)*size(vCin,2)) ~= 1
        error('ERROR: !!! Input argument into combineDims() is in wrong format - not a cell vector.')
    end
    %%% %%%. % There is more to check for - for example that bin centers occur at binSize
    %%% intervals, but i just assume that for now.
    
    nVs = length(vCin);
    
    %%% sum the lengths of the input vectors
    lengthOFvC_iul = 0;
    for v = 1:nVs
        if size(vCin{v},1)+size(vCin{v},2)-(size(vCin{v},1)*size(vCin{v},2)) ~= 1
            error('ERROR: !!! input bin centers not as vectors?')
        end
        lengthOFvC_iul = lengthOFvC_iul + length(vCin{v});
    end
    %%% %%%.
    
    %%% Initialise vC at initial upper limit (iul) size.
    vC = zeros(nVs + 1,lengthOFvC_iul);
    %%% %%%.
    
    %%% Fill in vC with all the values in the inputs, keeping a record of
    %%% from which input each value came.
    curLoc = 0;
    for v = 1:nVs
        vC((nVs + 1),(curLoc+1):(curLoc + length(vCin{v}))) = vCin{v};
        vC(v,(curLoc+1):(curLoc + length(vCin{v}))) = 1;
        curLoc = curLoc + length(vCin{v});
    end
    %%% %%%.
    
    %%% Sort the values in vC
    [~,I] = sort(vC((nVs + 1),:));
    vC = vC(:,I);
    %%% %%%.
    
    
    vC(nVs+1,:) = correctRounding(vC(nVs+1,:));
    %%% Remove duplicate values in vC, and maintain the record of from which input each value came.
    for k = 2:lengthOFvC_iul
        if vC((nVs + 1),k) == vC((nVs + 1),k-1)
            vC(1:nVs,k) = vC(1:nVs,k-1) + vC(1:nVs,k);
            vC(:,k-1) = 0;
        end
    end
    
    

    zeroLocs = (sum(vC) == 0);
    vC = correctRounding(vC(:,zeroLocs == 0));
    %%% %%%.
    % DONE.
end