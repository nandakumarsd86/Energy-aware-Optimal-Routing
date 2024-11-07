function[sol] = check_obj(max,sol)
B = 1:max;
a = round(sol);
A = unique(a);
if length(A) ~= length(a)    
   C = setdiff(B,A);
    m = 1;
    for k = 1:length(a)
        n = find(a == a(k));
        if length(n) ~= 1
            for j = 2:length(n)
                a(n(j)) = C(m);
                m = m+1;
            end        
        end
    end
end
sol = a;
end