function rnd = rnd_negbinom(mu, alpha, m, n)

r = 1/alpha;
p = r/(r+mu);

rnd = zeros(m,n);
for i=1:m
    for j = 1:n
        rnd(i,j) = negative_binomial(r, p);
    end
end

end

function rnd = negative_binomial(r, p)

shape = r;
scale = (1-p)/p;
Y = scale .* randg(shape);

rnd = poissrnd(Y);
end
