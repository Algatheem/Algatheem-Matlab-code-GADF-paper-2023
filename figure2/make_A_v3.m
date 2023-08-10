function A = make_A_v3(k, ell, nu, eta, beta, alpha, B, U, gamma, N )

M  = 4*N+2 ;
Mm = M - 1 ;
A  = zeros(M, M) ;

Bver    = B    * cos(gamma) ;
Bhor    = B    * sin(gamma) ;
betaver = beta * cos(alpha) ;
betahor = beta * sin(alpha) ;

% complete advection / diffusion terms for vorticity, then vector potential

n = - N ;
for i = 1:2:Mm
    fac0 = (n  +ell)^2 + k^2 ;
    facm = (n-1+ell)^2 + k^2 ;
    facp = (n+1+ell)^2 + k^2 ;
    for j = 1:M
        if j == i
            A(i,j) = - nu * fac0 - 1i  * (n+ell) * U;
        elseif j == i-2
            A(i,j) = 0.5 * k * ( - 1 + 1 / facm ) ;
        elseif j == i+2
            A(i,j) = 0.5 * k * (   1 - 1 / facp ) ;
        end
    end
    n = n + 1 ;
end

n = - N ;
for i = 2:2:M
    fac0 = (n  +ell)^2 + k^2 ;
%   facm = (n-1+ell)^2 + k^2 ;
%   facp = (n+1+ell)^2 + k^2 ;
    for j = 1:M
        if j == i
            A(i,j) = - eta * fac0 - 1i  * (n+ell) * U; 
        elseif j == i-2
            A(i,j) = - 0.5 * k ;
        elseif j == i+2
            A(i,j) =   0.5 * k ;
        end
    end
    n = n + 1 ;
end

if beta > 0   % include beta term
    n = - N ;
    for i = 1:2:Mm
        fac0 = (n  +ell)^2 + k^2 ;
 %      facm = (n-1+ell)^2 + k^2 ;
 %      facp = (n+1+ell)^2 + k^2 ;
        for j = 1:M
            if j == i
                A(i,j) = A(i,j) + 1i * ( - k * betahor  + (n+ell) * betaver ) / fac0 ;
            end
        end
        n = n + 1 ;
    end
end

if abs(Bver) > 1e-10   % include vertical B coupling terms
    
    n = - N ;
    for i = 1:2:Mm
        fac0 = (n  +ell)^2 + k^2 ;
%       facm = (n-1+ell)^2 + k^2 ;
%       facp = (n+1+ell)^2 + k^2 ;
        for j = 1:M
            if j == i+1   % couple H_n into G_n
                A(i,j) =  A(i,j) + 1i * k * Bver * fac0 ;
            end
        end
        n = n + 1 ;
    end
    
    n = - N ;
    for i = 2:2:M
        fac0 = (n  +ell)^2 + k^2 ;
%       facm = (n-1+ell)^2 + k^2 ;
%       facp = (n+1+ell)^2 + k^2 ;
        for j = 1:M
            if j == i-1   % couple G_n into H_n
                A(i,j) =  A(i,j) + 1i * k * Bver / fac0 ;
            end
        end
        n = n + 1 ;
    end
end

if abs(Bhor) > 1e-10  % include horizontal field coupling terms
        
    n = - N ;
    for i = 1:2:Mm
        fac0 = (n  +ell)^2 + k^2 ;
        facm = (n-1+ell)^2 + k^2 ;
        facp = (n+1+ell)^2 + k^2 ;
        for j = 1:M
            if j == i+1   % couple H_n into G_n
                A(i,j) =  A(i,j) + Bhor * 1i * (n+ell) * fac0 ; 
            end
            if j == i-1   % couple H_{n-1} into G_n
                A(i,j) =  A(i,j) + Bhor * 0.5 * 1i * k / eta * (facm-1) ; 
            end
            if j == i+3   % couple H_{n+1} into G_n
                A(i,j) =  A(i,j) + Bhor * 0.5 * 1i * k / eta * (facp-1) ; 
            end
        end
        n = n + 1 ;
    end
    
    n = - N ;
    for i = 2:2:M
        fac0 = (n  +ell)^2 + k^2 ;
        facm = (n-1+ell)^2 + k^2 ;
        facp = (n+1+ell)^2 + k^2 ;
        for j = 1:M
            if j == i-1   % couple G_n into H_n
                A(i,j) =  A(i,j) + Bhor * 1i * (n+ell) / fac0 ;
            end
            if j == i-3   % couple G_{n-1} into H_n
                A(i,j) =  A(i,j) + Bhor * 0.5 * 1i * k / ( eta*facm ) ;
            end
            if j == i+1   % couple G_{n+1} into H_n
                A(i,j) =  A(i,j) + Bhor * 0.5 * 1i * k / ( eta*facp ) ; 
            end
        end
        n = n + 1 ;
    end
    
end

end

