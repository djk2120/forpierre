function  [A2,ci]=get_A(j,gsc_mol)
    %get_A, solves for photosynthesis a la Medlyn 2011 et al.
    %   returns
    %     A2, gross photosynthesis [umol/m2/s]
    %   inputs
    %     j, [umol/m2/s]
    %     gsc_mol, stomatal conductance to CO2 [mol/m2/s]

    gam = 50;
    Ca = 400  ; %ppmv

    ci = 0;
    jump = 97;
    ct = 0;
    go = 1;

    while go
        ct = ct+ 1;
        ci = ci+ jump;

        A1 = j/4*(ci-gam)/(ci+2*gam);
        A2 = gsc_mol*(Ca-ci);

        if abs(A1-A2)<0.0001
            go   = 0;
        elseif A1>A2
            ci   = ci-jump;
            jump = jump/2.75;
        end

        if ct>100
            disp('photosynthesis not converging')
            go=0;
        end

    end

end

