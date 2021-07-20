

% File for flow laws storage

%table 21.1 gerya, from Ranalli 
MAT(2).pre_disc=2*2.5*1e4*((1e6)^-3.5);% dry olivine, Gerya table 21.1 Ranalli 1995
MAT(2).Edisc=532e3;
MAT(2).ndisc=3.5;

%issues with rising sides
%%
%table 1 from Bickert et al.(1) derived from Goetze 1978

MAT(2).pre_disc=2*7*1e4*((1e6)^-3);
MAT(2).Edisc=520e3;
MAT(2).ndisc=3;

%worse for sides rising
%%
%table 1 from Bickert et al, (2), derived from Karato 1986

MAT(2).pre_disc=2*1.10*1e5*((1e6)^-3.5);
MAT(2).Edisc=530e3;
MAT(2).ndisc=3.5;

%same issues as Ranalli not much of a difference
%%
%from Hirth and kohlstedt 2003, wet olivine

MAT(2).pre_disc=(2^3.5)*3.77e-14;%hirth and kohlstedt 2003
MAT(2).Edisc=520e3;
MAT(2).ndisc=3.5;
