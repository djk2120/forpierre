function vpd = get_vpd( rh,T )
% T [K], rh [0,1], vpd [Pa]

    lv  = 2.5e6  ;
    rv  = 461    ;
    es  = 611*exp(lv/rv*(1/273.16-1/T));
    vpd = (1-rh)*es;



end

