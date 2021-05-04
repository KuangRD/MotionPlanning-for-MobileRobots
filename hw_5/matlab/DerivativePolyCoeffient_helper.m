function poly_coefficients_array = DerivativePolyCoeffient_helper(n_order, rank, t)
    poly_coefficients_array = [];
    for i = 0:n_order
        if i >= rank 
            elem = t^(i-rank) * factorial(i)/factorial(i-rank);
        else
            elem = 0;
        end
        poly_coefficients_array = [poly_coefficients_array,elem];
    end
    
    %size(poly_coeffient_array);
end