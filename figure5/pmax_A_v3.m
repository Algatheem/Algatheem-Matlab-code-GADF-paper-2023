function [premax, pimmax, kmax] = pmax_A_v3(ell, nu, eta, beta, alpha, B, U, gamma, N )

lenk  = 200 ;
%spaced_out = linspace(-3,0,lenk) ;
%kvals   = 0.6 * 10.^spaced_out         ;

kvals = linspace(0.001,1.0,lenk) ;
prealvals = zeros(1,lenk) ; 
pimagvals = zeros(1,lenk) ; 

for j = 1:length(kvals)
    
    k = kvals(j) ; 
    
    A = make_A_v3(k, ell, nu, eta, beta, alpha, B, U, gamma, N ) ;
    
    p=eig(A) ;
    preal = real(p) ;
    pimag = imag(p) ;
    pmax=max(preal) ;
    
    k = find ( abs(preal - pmax) < 1e-10) ;
    
    if length(k) == 1
        k1 = k ;
    elseif imag(p(k(1))) > 0
        k1 = k(1) ;
    else
        k1 = k(2) ;
    end
    
    prealvals(j) = real(p(k1)) ;
    pimagvals(j) = imag(p(k1)) ; 

end

pmax = max(prealvals) ; 
k = find ( abs(prealvals - pmax) < 1e-10) ;
k = min(k) ; 

premax = prealvals(k) ; 
pimmax = abs(pimagvals(k)) ; 
kmax   = kvals(k) ; 
   
end