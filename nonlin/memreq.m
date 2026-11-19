
%	memory requirements for nonlin info structure

	a = 1;
	m = 6 * a;
	v = 4;
	n = 128;
	c = 128;

	r = [	n * m;
		n * v * c;
		n * m * v * c;
		n * m;
		n * v * c;
		n * (4*m*m+m) ]

	sum(r)
	100 * r / sum(r)

	w = 12;

	r = [	n * m;
		n * v * (1 + w);
		n * m * v * (1 + w);
		n * m;
		n * v * (1 + w);
		n * (4*m*m+m) ]

	sum(r)
	100 * r / sum(r)
