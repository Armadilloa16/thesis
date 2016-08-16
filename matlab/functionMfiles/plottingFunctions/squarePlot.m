function squarePlot(vx,vy)

    x_ill = min(vx);
    x_iul = max(vx);
    y_ill = min(vy);
    y_iul = max(vy);
    m = 5 + max((x_iul-x_ill),(y_iul-y_ill))/2;
    x_ll = ((x_ill+x_iul)/2) - m;
    x_ul = ((x_ill+x_iul)/2) + m;
    y_ll = ((y_ill+y_iul)/2) - m;
    y_ul = ((y_ill+y_iul)/2) + m;

    axis([x_ll x_ul y_ll y_ul])
    set(gca,'Ydir','reverse')
    
end