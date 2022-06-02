function xref = referenceorigins(m, n, p)

xref = zeros(3, m*n*p);
for i = 1:p 
    for j = 1:n 
        for k = 1:m 
            t = k + m*(j-1) + m*n*(i-1);
            xref(1, t) = k-1;  
            xref(2, t) = j-1;              
            xref(3, t) = i-1;              
        end
    end
end

end
