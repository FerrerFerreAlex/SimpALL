function [C,dC] = simp_all_interpolation(E1,E0,nu1,nu0)
% SIMP-ALL interpolation function by Alex Ferrer, March 2019.
%
% E0, nu0: Young modulus and Poisson ratio of the isotropic weak material
% E1, nu1: Young modulus and Poisson ratio of the isotropic strong material
%
% C(rho): Isotropic constitutive tensor interpolated from the isotropic
% material properties of the strong (E1,nu1) and the weak material (E0,nu0)
% through the "density function" rho.
%
% dC(rho): Derivative of the isotropic material constitutive tensor respect
% to the density function rho.
%
% Remark 1: No conditions are required to E1,E0,nu1,nu0 (except E>0 and
% -1<nu < 0.5), so we could have referred to the strong and weak material as
% material A and B. However, we prefer to name them strong and weak to
% follow topology optimization naming.
%
% Remark 2: This version of the SIMP-ALL interpolation function is intended
% for 2D plane stress problems.
%
% Remark 3: SIMP-ALL interpolation evidenced to fulfill the Hashin-Shtrikman
% (HS) bounds for all possible combination of the material properties.
% Thus, when the density function takes intermediate values, the resulting
% constitutive tensor can be interpreted as an homogenized constitutive
% tensor of a micro-structure made of the strong and weak material with a
% relative fraction volume of the intermediate density value.
%
% Remark 4: SIMP-ALL interpolation is free of penalization and heuristic
% parameters.
%
% Remark 5: The interpolation is built on the one hand such that the
% constitutive tensor coincides with C0 (constitutive tensor with E0,nu0
% material properties) in rho=0 and with C1 (respectively E1,nu1) in rho=1
% and on the other hand such that the constitutive tensor derivative in
% rho=0 and rho=1 coincides with the topological derivative of the
% compliance dC0 and dC1, i.e.,
%                  C(0)= C0,        C(1)= C1,
%                 dC(0)= dC0,      dC(1)= dC1.
%
% Remark 6: We opt for providing dC since when using the compliance
% (integral f*u = integral e(u):C(rho):e(u)) in a  topology optimization
% problem, its gradient is directly obtained as -e(u):dC(rho):e(u).
%
% Shear and bulk modulus
mu = @(E,nu) E/(2*(1+nu));
kappa = @(E,nu) E/(2*(1-nu));
%
mu0 = mu(E0,nu0);
mu1 = mu(E1,nu1);
%
kappa0 = kappa(E0,nu0);
kappa1 = kappa(E1,nu1);
%
% Auxiliar material property
eta_mu = @(mu,kappa) (kappa*mu)/(2*mu+kappa);
eta_kappa = @(mu,kappa) mu;
%
eta_mu0 = eta_mu(mu0,kappa0);
eta_mu1 = eta_mu(mu1,kappa1);
%
eta_kappa0 = eta_kappa(mu0,kappa0);
eta_kappa1 = eta_kappa(mu1,kappa1);
%
% Isotropic fourth order tensor in Voigt notation
I1 = eye(3);I1(3,3) = 1/2;
I2 = [1 1 0; 1 1 0; 0 0 0];
Aiso = @(alpha,beta) alpha*I1 + beta*I2;
%
% Coeficients (n= numerator, d = denominator) of the rational function
n01 = @(f0,f1,eta0,eta1) -(f1 - f0)*(eta1 - eta0);
n0 = @(f0,f1,eta0) f0*(f1 + eta0);
n1 = @(f0,f1,eta1) f1*(f0 + eta1);
d0 = @(f1,eta0) (f1 + eta0);
d1 = @(f0,eta1) (f0 + eta1);
%
n01_mu    = n01(mu0,   mu1,   eta_mu0,   eta_mu1);
n01_kappa = n01(kappa0,kappa1,eta_kappa0,eta_kappa1);
%
n0_mu    = n0(mu0,   mu1,   eta_mu0);
n0_kappa = n0(kappa0,kappa1,eta_kappa0);
%
n1_mu    = n1(mu0,   mu1,   eta_mu1);
n1_kappa = n1(kappa0,kappa1,eta_kappa1);
%
d0_mu    = d0(mu1,   eta_mu0);
d0_kappa = d0(kappa1,eta_kappa0);
%
d1_mu    = d1(mu0,   eta_mu1);
d1_kappa = d1(kappa0,eta_kappa1);
%
% Density function (symbolic in order for further differentiation)
rho = sym('rho','positive');
%
% SIMP-ALL as a rational function
f = @(n01,n0,n1,d0,d1,rho) ...
     (n01*(1-rho)*(rho) + n0*(1-rho) + n1*rho)/(d0*(1-rho)+d1*rho);
%
mu =    f(n01_mu,   n0_mu,   n1_mu,   d0_mu,   d1_mu,   rho);
kappa = f(n01_kappa,n0_kappa,n1_kappa,d0_kappa,d1_kappa,rho);
%
% Isotropic constitutive tensor.
Csym = simplify(Aiso(2*mu,kappa - mu));
%
% Shear and bulk modulus derivatives
dmu = diff(mu);
dkappa = diff(kappa);
%
% Isotropic constitutive tensor derivative
dCsym = Aiso(2*dmu,dkappa-dmu);
%
% In order to obtain a handle function instead of a symbolic expression and
% to directly evaluate the constitutive tensor (and derivative) for any
% vector of densities, we make use of matlabFunction().
nstre = 3; % Number of stress components in Voigt notation in plane stress
dC = cell(nstre,nstre);
for i=1:nstre
    for j=1:nstre
            if Csym(i,j)==0
                C{i,j} =  @(rho) zeros(size(rho));
                dC{i,j} = @(rho) zeros(size(rho));
            else
                C{i,j} = matlabFunction(Csym(i,j));
                dC{i,j} = matlabFunction(dCsym(i,j));
            end
    end
end
end

