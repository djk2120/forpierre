
clear

%parameter setup1
kmax_other_units = 220;
kmax  = kmax_other_units/1e6*18;
kmax_kg          = kmax_other_units*18/1000/1000;
z     = 30;
zr    = 3;
p1    = -1;
p2    = -4;
p50   = -1.5;
a     = 6;
param = [kmax,z,p1,p2,p50,a];


%meteorological conditions (constant)
R     = 450;
j     = 145;
T     = 305;
ga    = -1;
p_atm    = 1e5;
r_gas    = 8.314;
x        = 1.6*r_gas*T/p_atm;


%meteorological conditions (varying)
pv    = 0:-0.4:-2;      %soil potential (MPa)
rv    = 0.5:0.01:0.99;  %relative humidity (-)
nx    = length(rv);
ny    = length(pv);

%output variables
vpd   = zeros(nx,1);
pleaf = zeros(nx,ny);
A     = zeros(nx,ny);


for i = 1:nx
    rh = rv(i);
    vpd(i) = get_vpd(rh,T);
    gsc_max = medlyn(get_vpd(rh,T)/1000,j);
    gsw_max = gsc_max*x;
    
    for ii = 1:ny
        psoil = pv(ii);
        pen = [R,rh,T,gsw_max,ga];
        [q,pleaf(i,ii),gsw] = get_vwp(psoil,param,pen);
        A(i,ii)             = get_A(j,gsw/x);
    end
end

subplot(1,1,1)
plot(vpd/1000,A,'Linewidth',1.5)
ylim([0,25])

