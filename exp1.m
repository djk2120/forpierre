close all
clear

%parameter setup1
kmax_other_units = 220;
kmax  = kmax_other_units/1e6*18;
kmax_kg          = kmax_other_units*18/1000/1000;
z     = 30;
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
ci    = zeros(nx,ny);


for i = 1:nx  %cycle through vpd values
    rh = rv(i);
    
    %solve for maximum stomatal conductance
    %irrespective of psi_soil
    vpd(i) = get_vpd(rh,T);
    gsc_max = medlyn(get_vpd(rh,T)/1000,j);
    gsw_max = gsc_max*x;
    
    for ii = 1:ny  %cycle through psoil values
        psoil = pv(ii);
        pen = [R,rh,T,gsw_max,ga];  %arrange the forcing data correctly
        [q,pleaf(i,ii),gsw] = get_vwp(psoil,param,pen); %solve for leaf potential
        [A(i,ii),ci(i,ii)]  = get_A(j,gsw/x);  %solve for attenutated GPP
    end
end



%plotting
ll = cell(ny,1);
for i = 1:ny
    if i==1
        ll{i} = ['\psi_s=',num2str(round(pv(i),1)),' MPa'];
    else
        ll{i} = [num2str(round(pv(i),1)),' MPa'];
    end
    
end

xdk = figure('units','inches','position',[2,2,8,3],...
    'PaperSize',[8,3]);

subplot('Position',[0.1,0.15,0.3,0.8])
plot(vpd/1000,A,'Linewidth',1.5)
ylim([0,25])
xlabel('VPD (kPa)')
ylabel('GPP (g/m2/d)')
xlim([0,2.5])


subplot('Position',[0.5,0.15,0.475,0.8])
plot(vpd/1000,ci,'Linewidth',1.5)
xlabel('VPD (kPa)')
ylabel('C_i (ppmv)')
ylim([0,400])
legend(ll,'Location','eastoutside')
xlim([0,2.5])

print(xdk,'fig1','-dpdf')


