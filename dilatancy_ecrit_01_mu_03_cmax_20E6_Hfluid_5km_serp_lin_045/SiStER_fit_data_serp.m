%% determination of a reaction rate law for serpentinization modelling
clear all
close all
%% data from Malvoisin et al. (2012) (only one grain size for which the data are the best
Th = [350 300 270 250] + 273.15; % temperature in K
Ph = 500*ones(1,4); % pressure in bars
th = [6000 6020 5960 5180].*3600; % duration 
xih = [0.1 0.345 0.32 0.26]; % reaction progress measured
kh = -log(1-xih)./(th); % constant part of the reaction rate (dxi/dt = kh*(1-xi)) in s-1
Gs = [44 44 44 44]; % grain size (diameter) in microns

SSAh = 10.^(5.2-log10(Gs)); % from Brantley and Mellott: SSA0
SSA_ref = 10.^(5.2-log10(100)); % from Brantley and Mellott: to normalize the data acquired with Gs for a grain size of 100 microns
kh = kh./SSAh*SSA_ref; % kh for 100 microns

% other technique like the other
A =textread('dG_fo_sp_sup92_redlichkwong.tab'); % load dG data (2D)
Tfs = A(:,1); Tfs = reshape(Tfs,401,301);
Pfs = A(:,2); Pfs = reshape(Pfs,401,301);
dGfs = A(:,4); dGfs = reshape(dGfs,401,301)/2;
for i=1:size(Tfs,2) % fit dG data with a polynom (simpler for modelling afterwards)
    [a,b] = polyfit(Tfs(:,1),dGfs(:,i),1);
    coef_T(i) = a(1);
    clear a b
end
% final parameters of the fit of drG
alpha = mean(coef_T);
Tref = 367.8962985+273.15; % temperature of the equilibrium


% hydration only here
DGh =  (1-exp(((Th-Tref)*alpha)./(8.314*Th))); % for the hydration reaction
lambda_v = [1 2 3]; % exponent of the DrG function, TST, from Lasaga 1994 for 2 and Lasaga 1986 for the rest : 

for i = 1:length(lambda_v)
    DGst = (DGh).^lambda_v(i);
    [coeff1 rtemp]= fit(1./(8.314.*Th'),log(kh'./DGst'),'poly1')
    Ea_v(i) = -coeff1.p1
    k_v(i)  = exp(coeff1.p2);
    rs (i) = rtemp.rsquare;
end

Ea= Ea_v(rs == max(rs));
k = k_v(rs == max(rs));
lambda  = lambda_v(rs == max(rs));

% plot fit and data
figure
plot(Th,kh,'k*')
hold on
T_l = [5:1:450]+273.15;
k_c = k*exp(-Ea./(8.314.*T_l)).*((1-exp(((T_l-Tref)*alpha)./(8.314.*T_l)))).^lambda;
k_c(T_l > Tref) = 0; 
plot(T_l,k_c)
xlabel('Temperature (K)')
ylabel('k (s^{-1})')


