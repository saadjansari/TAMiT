function x_hat = Modulate(n0, x_hat)

[N, D] = GetDimensions(x_hat);
for d = 1:D
	x_hat = FilterSep(x_hat, d, exp(2i*pi*n0(d)*(0:N(d)-1)/N(d)));
end