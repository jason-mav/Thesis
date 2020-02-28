function [D, B] = matrix_balance(A)
    
% Numerical Recipes p.593 - Balancing (vedi anche articolo James)
% Trasformazione
% D^-1 * A * D = B 
    
    n = size(A,1);
    B = A;
    D = eye(n);
    
    RADIX = 2;
    sqrdx = RADIX*RADIX;
    
    done = 0;
    
    while done ~= 1 
        
        done = 1;
        
        for i = 1 : n 
            
            r = 0;
            c = 0;
            
            for j = 1 : n
                if j ~= i 
                    c = c + abs(B(j,i));
                    r = r + abs(B(i,j));
                end
            end
                
            if c ~= 0 && r ~= 0 
                g = r/RADIX;
                f = 1;
                s = c + r;
                
                while (c<g) 
                    f = f*RADIX;
                    c = c*sqrdx;
                end
                
                g = r*RADIX;
                
                while c > g 
                    f = f/RADIX;
                    c = c/sqrdx;
                end
                
                if (c + r)/f < 0.95*s 
                    done = 0;
                    g = 1/f;
                    D(i,i) = D(i,i)*f;

                    for j = 1 : n 
                        B(i,j) = B(i,j)/f; 
                    end
                    
                    for j = 1 : n
                        B(j,i) = B(j,i)*f;
                    end
                end
            end
        end
    end

end


