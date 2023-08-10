
clear all
close all

N=16;
rho=1;
mu=1;
U=0.0;
nu=0.4;
eta=nu;
epsi = nu ;
gamma = 0.0 ;
beta = 0.0 ;
alpha = 0.0 ;
kvals=linspace(0.001,0.9,200);
ell = 0.0 ;

Bs   = [ 0.0, 0.05, 0.10, 0.15, 0.20, 0.25 ]  ;

for j = 1:6

    B   = Bs(j) ;

    if j == 1
        col = '#0072BD'; %blue
    elseif j == 2
        col = '#A2142F'; %red
    elseif j == 3
        col = '#77AC30';%green
    elseif j == 4
        col = '#7E2F8E';%purple
    elseif j == 5
        col = '#EDB120';%orange
    elseif j ==6
        col = '#D95319'; %dark orange

    end
    for i=1:200

        k=kvals(i);

        A= make_A_v3(k, ell, nu, eta, beta, alpha, B, U, gamma, N ) ;

        p=eig(A) ;
        preal=real(p);
        pimag=imag(p);

        alfvenpreal = -k^2*(nu+eta)/2;
        alfvenpimag =1i*k*sqrt(B^2 - 0.25*k^2*(nu-eta)^2);
        %alfvenpimag =(1*1i/2)*sqrt((-k^4*((epsi/rho)+eta)^2)+(4*((epsi*eta/rho)*k^4+B^2*k^2/(rho*mu))));

        alfvenpvals(i) = alfvenpreal + alfvenpimag ;

        pmax=max(preal) ;


        k = find ( abs(preal - pmax) < 1e-10) ;

        if length(k) == 1
            k1 = k ;
        elseif imag(p(k(1))) > 0
            k1 = k(1) ;
        else
            k1 = k(2) ;
        end

        pmax = (p(k1)) ;

        pvals(i)=pmax;


    end

    % plotting

    figure(1)
    plot(kvals,real(pvals),'-','Color',col,'LineWidth',1.5)
    hold on
    %hold on
    hold on
    xlabel('$k$','Interpreter','LaTex')
    ylabel('$\mathrm{Re}\; p$','Interpreter','LaTex')
    %grid on
    xlim([0.0, 0.9])
    ylim([-0.04,0.06])
    set( gca                       , ...
        'FontName'   , 'Helvetica' , ...
        'FontSize'   , 20        );
    hold on



    figure(2)
    if j > 1
        plot(kvals, imag(pvals),'-','Color',col,'LineWidth',1.5)
        hold on
        plot(kvals,imag(alfvenpvals),'--','Color',col,'LineWidth',1.5)
        xlabel('$k$','Interpreter','LaTex')
        ylabel('$\mathrm{Im}\; p$','Interpreter','LaTex')
        xlim([0.0, 0.9 ])
        ylim([0.0, 0.03])    %grid on
        set( gca                       , ...
            'FontName'   , 'Helvetica' , ...
            'FontSize'   , 20        );
        hold on
    end
end

%%

figure(1)
plot(kvals,real(alfvenpvals),'k--','linewidth',1.5)

legend('$B_0=0$',   '$B_0=0.05$','$B_0=0.1$', '$B_0=0.15$', ...
           '$B_0=0.2$', '$B_0=0.25$', 'Alfven', ...
           'Interpreter','LaTex','Location','NorthEast')

figure(2)
plot(-kvals,-real(alfvenpvals),'k--','linewidth',1.5) % Dummy line
legend('$B_0=0.05$', '', '$B_0=0.1$', '', '$B_0=0.15$','',  ...
           '$B_0=0.2$', '', '$B_0=0.25$', '', 'Alfven', ...
           'Interpreter','LaTex','Location','NorthEast')

