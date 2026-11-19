% dbvs.m
% Matlab file for debugging vector spline smoothing algorithm
% and its mex file vs_mex()

% SETUP DATA
if (exist('ti') == 0)
	equispace = 0;
	ti = [1 3 4 8 16 17 19]';
	N = length(ti);
	if (equispace)
		ti = [1:N]';
	end
end

if (exist('yr') == 0)
	y1 = 10*[1 2 3 3  3  2  1]';
	yr = [y1, 10*[1:length(ti)]'];
	alpha = [-1 1]';
	M = length(alpha);
end

if (exist('yi') == 0)
	sig = [3 2; 2 4]/6;
	isig = inv(sig);
	Sig = sig(:) * ones(ti');
	iSig = isig(:) * ones(ti');
	Pi = kron_mex(eye(N), sig);

%	These next few lines make it non-iid
	iid = 0;
	if (iid == 0)
		sig2 = [0.01 0; 0 0.01];
		isig2 = inv(sig2);
		for n=[3,N]
			Pi((n-1)*M+[1:M],(n-1)*M+[1:M]) = sig2;
			Sig(:,n) = sig2(:);
			iSig(:,n) = isig2(:);
		end
	end
	iPi = inv(Pi);

	rand('seed', 0)
	Pi2 = sqrtm(Pi);
	yi = yr';
	yi(:) = Pi2 * rand(yi(:));
	yi = yr + yi';
end

%	diary data2, disp([N M]), disp(alpha'), disp(ti'), disp(yi), disp(iSig'), diary off

%
%	Mex smoothed version
%
	[Vur,Vcv,Vgcv,Vys] = vs_mex(alpha,ti,yi',iSig);

%
%	FIND SMOOTHED VERSION in MATLAB - ALA MEX FILE
%
	eyeMN = eye(M*N);
	[qt, r] = spline_mat(ti);
	Qt = kron_mex(qt,eye(M));
	rb = mat2band(r, 1);
	qtb = mat2band([qt; zeros(2, ncol(qt))], 2);
%	qtq = qt * qt';
%	qtqb = mat2band(qtq, 2);
%	qtq_sig = kron_mex(qtq, sig);
%	qtsqb = mat2band(qtq_sig, 3*M-1);
	qtsq = Qt * Pi * Qt';
	qtsqb = mat2band(qtsq, 3*M-1);

	yit = yi';
	QtY = Qt * yit(:);

	alpha_new = exp(alpha) * (N-1)/(ti(N)-ti(1));
	ra = kron_mex(r, diag(1 ./ alpha_new));
	rab = mat2band(ra, 2*M-1);
	b = ra + qtsq;
	bb = mat2band(b, 3*M-1);
	[l, d] = mychol(b);
	d = diag(d)';
	lb = mat2band(l', 3*M-1);

	ca = b \ QtY;
	c = zeros(M, N-2);
	c(:) = ca;
	c = diag(1 ./ alpha_new) * c;
	q_c = Qt' * ca;
	ys = zeros(M, N);
	ys(:) = yit(:) - Pi * q_c;
	ys1 = ys';

	ib = inv(b);
	ibb = mat2band(ib, 3*M-1);
	hat = eyeMN - Pi * Qt' * ib * Qt;

%	mask = kron_mex(eye(N), ones(M,M));
	f = kron_mex(qt',eye(M)) * ib * kron_mex(qt, eye(M));
	f = f .* kron_mex(eye(N), ones(M,M));
	[fl, fd] = mychol(f);
	fd = diag(fd);
	fi = inv(f);
	fipi = fi * iPi;
	err = yi' - ys1';
	werr = err;
	werr(:) = fipi * err(:);
	nerr = err; nerr(:) = iPi * werr(:);
	cv = sum(werr .* nerr);
	cv = cv / N;

	rss = err(:)' * iPi * err(:);
	errn = err(:)' * err(:);
	ur = (errn - 2 * trace(Pi * (eyeMN - hat)) + trace(Pi)) / N;
	gcv = N * rss / trace(eyeMN - hat)^2;
	if (equispace)
		roughness = c(:)' * kron_mex(r, diag(alpha_new)) * c(:)
	else
		roughness = -1;
	end

%	PENALTY
	penin = QtY' * (ra \ QtY)
	z = Qt * ys(:); penout = z' * (ra \ z)
