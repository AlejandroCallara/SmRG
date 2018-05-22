function rnd = rnd_negbinom(mu, alpha, m, n)
%|====================================================================================
%|RND_NEGBINOM Generate a MxN matrix with samples drawn from NEGBIN_MU distribution
%|
%|  RND = RND_NEGBINOM(MU,ALPHA,M,N) returns an array of random numbers chosen 
%|                                   from the Negative Binomial distribution with 
%|                                   parameters MU and ALPHA. 
%|                                   The size of RND is [M,N], with M,N scalars.
%| 
%|  See also nbinrnd_mu, nbinfit_mu, nbinpdf_mu, nbinlike_mu, random.
%|
%|  Last revision:
%|  22 May 2018
%|  Michele Scipioni, University of Pisa
%|
%|====================================================================================


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
