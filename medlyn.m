function gs= medlyn(vpd,j)
    %medlyn, solves for medlyn stomatal conductance
    %   returns
    %    gs, stomatal conductance [mol/ms/s]
    %   inputs
    %    vpd, [kPa]
    %    j,   [umol/m2/s]
    
    g0  = .001       ; % mol / m2 / s
    g1  = 5       ; % kPa ^ 0.5 

    Ca = 400  ; %ppmv

    go = 1;
    ct = 0;

    oldA  = 10   ; %umol/m2/s
    while go
        ct = ct+1;

        gs = g0 + (1+g1/vpd^(0.5))*oldA/Ca  ; %mol/m2/s
        A  = get_A(j,gs);

        if abs(A-oldA)<.001
            go=0;
        else
            oldA=A;
        end

        if ct>100
            go=0;
            disp('medlyn not converging')
        end
    end
            


end

