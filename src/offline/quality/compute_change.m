function d = compute_change(a, b)

EPSEPS = 1.0e-3;

l = numel(a);
d = zeros(size(a));
for i = 1 : l
    d(i) = abs(a(i) - b(i)) / max([abs(a(i)), abs(b(i)), EPSEPS]);
end

end