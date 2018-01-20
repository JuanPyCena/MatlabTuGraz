function A = getA(w_plus_delta,w_0, delta_x)
    rows = length(w_0);
    cols = length(delta_x);
    
    A = zeros(rows, cols);
    
    for i = 1:rows
       %w component
       for j = 1:cols
          % derivative by parameters
          if delta_x(j) == 0
              delta_x(j) = 1;
          end
          A(i,j) = (w_plus_delta(i,j) - w_0(i)) / delta_x(j);
          
%           if A(i,j) < minimum
%               A(i,j) = 0;
%           end    
       end
    end
        
end

