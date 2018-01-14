function [vec] = explicit(vec,dt)
A_3 = 5e3;
db = 1000;
t_krit = db^2/(2*A_3);
t_max = 172800;
N_anf = 2e10;
N_end = 20e10;
index_max = 251;

for t = 1:t_max /(t_krit)
    for j = 2:size(vec,1)-1
        delta_N(j-1,t) = (dt * t_krit * A_3 / db^2) * (vec(j+1,t) - 2*vec(j,t) + vec(j-1,t));
    end
    
    vec(1,t+1) = N_anf;
    vec(size(vec,1),t+1) = N_end;
    vec(2:end-1,t+1) = vec(2:end-1,t) + delta_N(:,t);
    
    %Nicht unter Randwerte fallen lassen
    i_max = 251;
    index_a = find(vec(1:i_max,t+1) < N_anf);
    vec(index_a,t+1) = N_anf;
    index_e = find(vec(i_max:end,t+1) < N_end);
    vec(index_e + index_max-1,t) = N_end;
end
end

