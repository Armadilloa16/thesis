function fnam = fnamSelect(folnam,X,Y)

    fnam = [folnam '/peaklists/0_R00X'];

    if X > 99
        fnam = [fnam num2str(X) 'Y'];
    elseif X > 9
        fnam = [fnam '0' num2str(X) 'Y'];
    else
        fnam = [fnam '00' num2str(X) 'Y'];
    end

    if Y > 99
        fnam = [fnam num2str(Y)];
    elseif Y > 9
        fnam = [fnam '0' num2str(Y)];
    else
        fnam = [fnam '00' num2str(Y)];
    end
    
    fnam = [fnam '.txt'];
    

end