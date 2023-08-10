clear all
close all

% does contour plot (eps, B0) 
 
ell     = 0.0 ; 
beta    = 0.0 ; 
alpha   = 0.0 ; 
gamma   = 0.0 ;
Prandtl = 0.5 ; 
U       = 0.0 ; 
 
for N = 8
    
% do an eps - B0 growth rate plot

nuvals  = linspace(0.01, 0.4, 100) ;
B0vals  = linspace(0.3,  0.8, 100) ;

for i = 1:length(nuvals)
    disp([i])
    for j = 1:length(B0vals)
        
        nu  = nuvals(i) ; 
        eta = 1/Prandtl*nuvals(i) ; 
        B   = B0vals(j) ; 
        
  %      pmaxvals(j,i) = pmax_A_fminbd_v1(ell, nu, eta, beta, alpha, B, gamma, N ) ;
  
  %      [ premaxvals(j,i), pimmaxvals(j,i), kmaxvals(j,i)] ...
  %            = pmax_A_theory_1_v2(ell, nu, eta, beta, alpha, B, gamma, N ) ;
  
        [ premaxvals(j,i), pimmaxvals(j,i), kmaxvals(j,i)] ...
                      = pmax_A_v3(ell, nu, eta, beta, alpha, B, U, gamma, N ) ;

    end
end

% save('run_1_data')

%% post process for contour maps

%load('run_1_data')

for i = 1:length(nuvals)
    for j = 1:length(B0vals)
        if premaxvals(j,i) < 0
            premaxvals(j,i) = 0.0  ;
        end
    end
end

maxcolmap = max(max(premaxvals))/255 ;

for i = 1:length(nuvals)
    for j = 1:length(B0vals)
        if premaxvals(j,i) > 0 && premaxvals(j,i) < maxcolmap
            premaxvals(j,i) = maxcolmap  ;
        end
    end
end

nustar = sqrt(Prandtl/2  * (1-Prandtl)/(1+Prandtl)) 

%%
close all

figure(10*N/4+4)
cmap = colormap('jet') ;
cmap(1,1:3) = 0.0 ;
colormap(cmap)
pcolor(nuvals, B0vals, premaxvals) 
shading flat
hold on 
cvals = [ 1e-19, 1e-19 ] ;
contour(nuvals, B0vals, premaxvals , cvals, 'w', 'LineWidth', 1.2)
plot([nustar, nustar], [0.0 10.0], 'w:', 'LineWidth',1.2)
xlabel('$\nu$','Interpreter','LaTex')
ylabel('$B_0$','Interpreter','LaTex')
title('$\mathrm{Re}\; p_{\max}$','Interpreter','LaTex')
colorbar
set( gca                       , ...
    'FontName'   , 'Helvetica' , ...
    'FontSize'   , 20          );
shading flat

figure(10*N/4+5)
shading flat
cmap = colormap('jet') ;
cmap(1,1:3) = 0.0 ;
colormap(cmap)
pcolor(nuvals, B0vals, pimmaxvals) 
shading flat
hold on 
%contour(nuvals, B0vals, premaxvals , cvals, 'w')
xlabel('$\nu$','Interpreter','LaTex')
ylabel('$B_0$','Interpreter','LaTex')
title('$\mathrm{Im}\; p_{\max}$','Interpreter','LaTex')
colorbar
set( gca                       , ...
    'FontName'   , 'Helvetica' , ...
    'FontSize'   , 20          );

figure(10*N/4+3)
shading flat
colormap('jet')
pcolor(nuvals, B0vals, kmaxvals) 
shading flat
hold on 
contour(nuvals, B0vals, kmaxvals , cvals, 'w')
xlabel('$\nu$','Interpreter','LaTex')
ylabel('$B_0$','Interpreter','LaTex')
title('$k_{\max}$','Interpreter','LaTex')
colorbar
set( gca                       , ...
    'FontName'   , 'Helvetica' , ...
    'FontSize'   , 20          );


drawnow
hold on 

end

 