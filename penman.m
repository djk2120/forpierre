function q = penman(R,rh,T,gs,ga)
% penman, a function to calculate penman ET
    % R  = W/2
    % rh = [0,1]
    % T  = K
    % gs = m/s
    % ga = m/s, default = 1/20

    if ga<0
        ga = 1/20;
    end
    
    vpd = get_vpd(rh,T);
    es  = vpd/(1-rh);
    
    lv  = 2.5e6  ; % J/kg
    rv  = 461    ; % J/kg/K


    P   = 1e5    ; % Pa
    eps = 0.622  ;
    rho = (P-(es*rh))/(eps*rv*T)+(es*rh)/(rv*T);

    cp  = 1004   ; % J/kg/K
    d   = lv/rv*es/T^2;
    g   = P*cp/(eps*lv);

    dq = eps*vpd/P;

    le = (d/g*R+rho*lv*ga*dq)/(1+d/g+ga/gs);
    q  = le/lv; %mm/s


end

