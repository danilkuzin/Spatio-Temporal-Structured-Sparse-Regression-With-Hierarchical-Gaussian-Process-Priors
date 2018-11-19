function r = max_change(a, b)

EPSEPS = 1.0e-3;

[l, m] = size(a);

r = 0;

for i = 1 :l
   for j = 1 : m
       s = abs(a(i, j) - b(i, j)) / max([abs(a(i, j)), abs(b(i, j)), EPSEPS]);
       if (s > r)
           r = s;
       end
   end
end

end
