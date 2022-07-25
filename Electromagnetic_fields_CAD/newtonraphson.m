function [n_iterazioni, zero] = newtonraphson(f, dfdx, nitermax, x0, tol, flag)

%Metodo iterativo per il calcolo dello zero di una funzione non lineare f.

ftol = abs(f(x0) * tol);
xtol = abs(f(x0) * tol);
xk = x0;
fk = f(xk);

for k = 1:nitermax
    
    switch (flag)
        case 1
            derivata = dfdx(xk);
        case 2
            delta = 0.00001;
            valass1 = abs(delta * xk);
            valass2 = abs(xk);
            den = 2 * delta * valass2;
            derivata = (f(xk + valass1) - f(xk - valass1))/den;
    end
    
    dxk = -fk/derivata;
    xk = xk + dxk;
    
    if (fk < ftol) && (dxk < xtol)
        n_iterazioni = k;
        zero = xk;
        return
    end
    
    fk = f(xk);
    
end

if (fk > ftol) || (dxk > xtol)
        error('Non è stata raggiunta la convergenza entro il numero di iterazioni massimo inserito')
end
    
end