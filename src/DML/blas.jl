
function densematmul(A, B)

    mA = size(A,1);
    nA = size(A,2);
    mB = size(B,1);
    nB = size(B,2);

    if nA != mB
        error("Input matrices are not multiplable")
    end

    C = zeros(mA,nB);
    for i = 1:mA
        for j = 1:nB
            for k = 1:nA
                C[i,j] = C[i,j] + A[i,k]*B[k,j];
            end
        end
    end

    return C

end

function densesparsematmul(A, B)

    mA = size(A,1);
    nA = size(A,2);
    mB = size(B,1);
    nB = size(B,2);

    if nA != mB
        error("Input matrices are not multiplable")
    end

    C = zeros(mA,nB);
    for j = 1:nB
        for k = 1:mB
            if abs(B[k,j]) > 0
                for i = 1:mA    
                    C[i,j] = C[i,j] + A[i,k]*B[k,j];
                end
            end
        end
    end

    return C
end

function densesparsematmul(V, W, A)

    M = size(W,1);
    N = size(A,1);

    D = zeros(M,N);
    for n = 1:N    
        for j = 1:N    
            if abs(A[n,j]) > 0
                for m = 1:M        
                    for i = 1:M
                        D[m,n] = D[m,n] + V[i,j]*W[i,m]*A[n,j];
                    end
                end
            end
        end
    end
    
    return D
end


    
    
    
    
    
    
    
    