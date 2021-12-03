function s_new = argo_sbe_CPcor (p,t,s,cpcor_new)
% function s_new = argo_sbe_CPcor (p,t,s,cpcor_new)
% Calculat new salinity using new CPcor (cpc_new)
% input: vector p,t,s, coeff. cpc
% output: vector s_new  (salinity)

% Coefficients
ctcor     = 3.25e-6;
cpcor_old = -9.57e-8;

% Calculate conductivity
crat = sw_cndr (s,t,p);
c = crat.*sw_c3515;

% Calculate new conductivity
n_old = 1.0 + ctcor.*t + cpcor_old.*p;
n_new = 1.0 + ctcor.*t + cpcor_new.*p;
c_new = c.*n_old./n_new;

% Calculate new salinity
s_new = sw_salt (c_new./sw_c3515,t,p);
end

