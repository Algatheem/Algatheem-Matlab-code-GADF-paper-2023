clear all
close all 

%  do a contour plot of some case or other! 

k     = 0.4 ; 
ell   = 0.0 ; 
nu    = 0.4 ; 
eta   = nu ; 
beta  = 0.0 ; 
alpha = 0.0 ; 
B     = 0.25 ; 
U     = 0.0 ; 
gamma = 0.0 ; 
N     = 10 ; 

A = make_A_v3(k, ell, nu, eta, beta, alpha, B, U, gamma, N ) ; 

p=eig(A) ;
[V, D] = eig(A) ; 

preal = real(p) ;
pimag = imag(p) ;
pmax=max(preal) ;

kind = find ( abs(preal - pmax) < 1e-10) ;

if length(kind) == 1
    kind1 = kind ;
elseif imag(p(kind(1))) > 0
    kind1 = kind(1) ;
else
    kind1 = kind(2) ;
end

preal = real(p(kind1)) ;
pimag = imag(p(kind1)) ;

pvec = V(:, kind1) ; 
Gn   = pvec(1:2:length(pvec)-1) ; 
Hn   = pvec(2:2:length(pvec)  ) ; 

% reconstruct the eigenfunction on a 100 x 100 grid 

xn = 101 ;
yn = 201 ; 
xxx = linspace(0,1,xn) * 2 * pi      ; 
yyy = linspace(0,1,yn) * 2 * pi / k  ;

F   = zeros(1,xn) ; 
G   = zeros(1,xn) ; 
H   = zeros(1,xn) ; 
J   = zeros(1,xn) ; 

[ XXX, YYY ] = meshgrid(xxx,yyy) ; 

% First reconstruct the functions of x 

for i = 1:xn
    for nnn = 1:2*N+1
        n    = nnn - 1 - N ; 
        lapfac = (n+ell)^2 + k^2 ; 
        
        F(i) = F(i) + Gn(nnn) * exp( 1i * ( n + ell) * xxx(i) ) / lapfac ; 
        G(i) = G(i) + Gn(nnn) * exp( 1i * ( n + ell) * xxx(i) ) ; 
        H(i) = H(i) + Hn(nnn) * exp( 1i * ( n + ell) * xxx(i) ) ; 
        J(i) = J(i) + Hn(nnn) * exp( 1i * ( n + ell) * xxx(i) ) * lapfac ; 

    end
end



% figure(13)
% 
% subplot(2,2,1)
% plot(xxx, real(F), 'b', xxx, imag(F), 'r', xxx, abs(F), 'g') 
% title('F') 
% subplot(2,2,2)
% plot(xxx, real(G), 'b', xxx, imag(G), 'r', xxx, abs(G), 'g') 
% title('G') 
% subplot(2,2,3)
% plot(xxx, real(H), 'b', xxx, imag(H), 'r', xxx, abs(H), 'g') 
% title('H') 
% subplot(2,2,4)
% plot(xxx, real(J), 'b', xxx, imag(J), 'r', xxx, abs(J), 'g') 
% title('J') 

% reconstruct full contour plot 

for i = 1:xn
    for j = 1:yn
        Ffield(i,j) = F(i) * exp(1i * k * yyy(j) )  ; 
        Gfield(i,j) = G(i) * exp(1i * k * yyy(j) )  ; 
        Hfield(i,j) = H(i) * exp(1i * k * yyy(j) )  ; 
        Jfield(i,j) = J(i) * exp(1i * k * yyy(j) )  ; 
    end
end

%% figures

figure(1)
colormap (parula)
shading flat 
contourf(XXX,YYY,real(Ffield)',15)
xlabel('$x$','Interpreter','LaTex')
ylabel('$y$','Interpreter','LaTex')
title('$\psi$','Interpreter','LaTex')
shading flat 
colorbar
axis on
xticks((0))
text(6,-0.8,'$2\pi$','interpreter','latex','fontsize',20)
yticks((0))
text(-0.6,15,'$\displaystyle{\frac{2\pi}{k}}$','interpreter','latex','fontsize',20)
set( gca                       , ...
    'FontName'   , 'Helvetica' , ...
    'FontSize'   , 20          );


figure(2)
colormap (parula)
shading flat 
contourf(XXX,YYY,real(Hfield)',15)
xlabel('$x$','Interpreter','LaTex')
ylabel('$y$','Interpreter','LaTex')
title('$a$','Interpreter','LaTex')
shading flat 
colorbar
axis on
xticks((0))
text(6,-0.8,'$2\pi$','interpreter','latex','fontsize',20)
yticks((0))
text(-0.6,15,'$\displaystyle{\frac{2\pi}{k}}$','interpreter','latex','fontsize',20)
set( gca                       , ...
    'FontName'   , 'Helvetica' , ...
    'FontSize'   , 20          );

% figure(14)
% 
% subplot(2,2,1)
% contour(xxx, yyy, real(Ffield)', 20)
% title('F - stream function') 
% subplot(2,2,2)
% contour(xxx, yyy, real(Gfield)', 20)
% title('G - vorticity') 
% subplot(2,2,3)
% contour(xxx, yyy, real(Hfield)', 20)
% title('H - vector potential') 
% subplot(2,2,4)
% contour(xxx, yyy, real(Jfield)', 20)
% title('J - current') 
% 
% 
% figure(15)
% 
% subplot(2,2,1)
% contour(xxx, yyy, imag(Ffield)', 20)
% title('F - stream function') 
% subplot(2,2,2)
% contour(xxx, yyy, imag(Gfield)', 20)
% title('G - vorticity') 
% subplot(2,2,3)
% contour(xxx, yyy, imag(Hfield)', 20)
% title('H - vector potential') 
% subplot(2,2,4)
% contour(xxx, yyy, imag(Jfield)', 20)
% title('J - current') 