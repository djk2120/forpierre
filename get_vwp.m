function [q1,pleaf,gsw] = get_vwp(psoil,param,pen)
%get_vwp, solves for equilibrium leaf water potential
%    returns
%     q1   ,  attenuated ET , [mm/s]
%     pleaf,  leaf potential, [MPa]
%     gsw  ,  stomatal conductance (to water), [m/s]

kmax  = param(1);
z     = param(2);
p1    = param(3);
p2    = param(4);
p50   = param(5);
a     = param(6);

R     = pen(1);
rh    = pen(2);
T     = pen(3);
gsw_max = pen(4);
ga    = pen(5);


grav  = z*9.8e-3;
step  = 0.5;
pleaf = psoil;
f   = (0.5*(psoil+pleaf)-p2)/(p1-p2);

if f<0.01
    q1 = 0;
    pleaf = psoil-grav;
    gsw = 0;
else

go = 1;
ct = 0;

while go
    pleaf = pleaf-step;
    ct = ct+1;
    
    f   = min(1,max(0,(0.5*(psoil+pleaf)-p2)/(p1-p2)));
    q1  = kmax/z*f*(psoil-pleaf-grav);
    
    h   = 1/(1+(pleaf/p50)^a);
    gsw = h*gsw_max;
    q2  = penman(R,rh,T,gsw,ga);
    
    if ct>50
        go = 0;
        if f<0.01
            'small f'
            q1 = 0;
            gsw = 0;
        else
            'vwp not converging'
            q1 = -1;
            
        end
    end
    
    if abs(q1-q2)<1e-9
        go = 0;
    elseif q1>q2||f==0
        pleaf = pleaf+step;
        step  = step/2.7;
    end
  
end
end

end
