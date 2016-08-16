function B = hcl2rgb(A)
    
        if length(size(A)) == 2
    
            H = A(:,1);
            C = A(:,2);
            L = A(:,3);

            Luv = A;
            Luv(:,1) = L;
            Luv(:,2) = cos((pi*H)/180).*C;
            Luv(:,3) = sin((pi*H)/180).*C;

            B = colorspace('Luv->RGB',Luv);
            
        elseif length(size(A)) == 3
            
            H = A(:,:,1);
            C = A(:,:,2);
            L = A(:,:,3);

            Luv = A;
            Luv(:,:,1) = L;
            Luv(:,:,2) = cos((pi*H)/180).*C;
            Luv(:,:,3) = sin((pi*H)/180).*C;

            B = colorspace('Luv->RGB',Luv);
            
        else
            error('ERROR: !!! Input to hcl2rgb is funky')
        end
    
end