function Use_example_of_simp_all_interpolation
% Use example of the SIMP-ALL interpolation function for computing the
%compliance and its gradient, by Alex Ferrer, May 2018.
%
% From a given material isotropic properties values
E1 = 1;
E0 = 0.01;
nu1 = 1/3;
nu0 = 1/3;
% Build the SIMP-ALL interpolation of the constitutive tensor as
[C,dC] = simp_all_interpolation(E1,E0,nu1,nu0);
%
%In the context of a FE problem, for a given number of gauss points ngaus
% and number of elements nelem of for example
nelem = 1000;
ngaus = 4;
%
% and for a given a density values 0 <= rho <= 1 of dimension
% dim(rho) = ngaus x nelem  of for example
rho = rand(ngaus,nelem);
%
% and for a given strain tensor values in Voigt notation (nstre = 3)
% obtained from thesolution of a finite element problem KU=F of for example
nstres = 3;
strain = rand(nstres,ngaus,nelem);
%
%Although it would be cheaper with F*U, defining dV (weight*Jacobian)
dV = rand(ngaus,nelem);
%the compliance can be computed as
c = 0;
for istres = 1:nstres
    for jstres = 1:nstres
      str = squeeze(strain(jstres,:,:));
      aux = str.*C{istres,jstres}(rho).*str.*dV;
      c = c + sum(aux(:));
    end
end
%
% and the gradient of the compliance as
g = zeros(ngaus,nelem);
for istres = 1:nstres
    for jstres = 1:nstres
        str = squeeze(strain(jstres,:,:));
        g = g - str.*dC{istres,jstres}(rho).*str;
    end
end
end
