function [CY, Y, TY, G, CvecG, M] = LocalPolyAnalv2(data, method);%%		Estimates the frequency response function (FRF) and the output noise covariance matrix from  %       arbitrary input/output data via a local polynomial least squares approximation of the %       plant transfer function and the plant and noise transient terms. The input signal may be%       exactly zero in (parts of) the frequency band of interest.%%       For nonlinear systems the FRF is the best linear approximation and the output covariance%       is the sum of the noise covariance and the covariance of the stochastic nonlinear distortions.  %%       If no input data is provided, then the algorithm simplifies to nonparametric %       time series analysis (noise power spectrum estimation).%%       Warning: the estimated frequency response matrix is meaningless in those %       frequency bands were the input is exactly zero.%       %%	function [CY, Y, TY, G, CvecG, M] = LocalPolyAnal(data, method);%%	Output parameters%%		CY      =	struct('n', [], 'm', [], 'm_nt', [])%                       CY.n	:   sample covariance matrix output noise, size ny x ny x F%                                       Usage: calculation uncertainty bounds on FRF%                       CY.m	:	sample covariance matrix sample mean output, size ny x ny x F%                                       Usage: weighting in frequency domain maximum likelihood estimation %                       CY.m_nt	:	sample covariance matrix sample mean output with transient removed, size ny x ny x F %                                       Usage: weighting in frequency domain maximum likelihood estimation %%		Y      =   struct('m', [], 'm_nt', [])%                       Y.m     :   sample mean output, size ny x F%                                       Usage: frequency domain maximum likelihood estimation %                       Y.m_nt  :   sample mean output with transient removed, size ny x F%                                       Usage: - frequency domain maximum likelihood estimation %                                              - leakage free output DFT spectra  %%		TY      =	sum plant and noise transient contribution at output, size ny x F %                       Usage: - transient removal in sample mean Ym(k) %                              - transient removal in output spectrum Y(k)%%       G       =   estimated frequency response matrix, size ny x nu x F%%       CvecG   =   covariance matrix vec(G), size (ny*nu) x (ny*nu) x F %%       M       =   equivalent number of independent experiments = number%                   of degrees of freedom in the residuals + 1 %%	Input parameters%%		data	=	structure containing the input/output data in the frequency band of interest %                   struct{'Y', [], 'U', [], 'f', [])%                       data.Y              =   output signal, size ny x F%                       data.U              =   input signal, size nu x F%                       data.freq           =   frequency vector in Hz or in DFT numbers, size 1 x F (optional)  %                                               default: [1:1:F]%%       method  =   structure containing the parameters of the method used (optional) %                   struct('order', [], 'moment', [], 'transient', [], 'step', [])%                       method.order        =	order of the polynomial approximation (optional; default 2) %                       method.moment       =	determines the existence of moment of the CY^-1 (optional; default 0)%                                                   moment = 0: CY is of full rank%                                                   moment = 1: then the expected value of CY^-1 exists%                                                   moment = 2: then the second order moments of CY^-1 exist %                       method.transient    = 	determines the estimation of the transient term (optional; default 1)  %                                                   1: transient term is estimates %                                                   0: no transient term is estimates %                       method.step         =   determines at which entries of data.freq the output parameters are calculated:   %                                               (optional; default 1)%                                                   data.freq(1:step:end) %% Rik Pintelon, July 2008% version 27 March 2009%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialisation variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%[ny, F] = size(data.Y);                                     % number of outputs ny, and number of frequencies F try    if isempty(data.freq)        data.freq = [1:1:F];    else        data.freq = data.freq(:).';        data.freq = data.freq/min(diff(data.freq));         % normalisation for improving the numerical conditioning      end % if   catch	data.freq = [1:1:F];    end % trytry    if isempty(method)        method = struct('order', 2, 'width', 1, 'moment', 0);    endcatch    method = struct('order', 2, 'width', 1, 'moment', 0);end % trytry	if isempty(method.order)		method.order = 2;	endcatch	method.order = 2;end % trytry	if isempty(method.width)		method.width = 1;	endcatch	method.width = 1;end % trytry	if isempty(method.moment)		method.moment = 0;	endcatch	method.moment = 0;end % trytry	if isempty(method.transient)		method.transient = 1;	endcatch	method.transient = 1;end % trytry    if isempty(method.step)        method.step = 1;    end % if   catch	method.step = 1;    end % trytry	if isempty(data.U)		data.U = [];	endcatch	data.U = [];endR = method.order;               % order polynomial methodwidth = method.width;           % frequency width in DFT samplestransient = method.transient;   % if 1 then transient is estimated (default); otherwise 0mm = method.moment;             % existence first or second order moment inv(CY) Fstep = method.step;            % frequency stepSelectFreq = [1:Fstep:F].';     % entries of data.freq at which the output parameters are calculated Fselect = length(SelectFreq);   % number of frequencies at which the output parameters are calculated  nu = size(data.U, 1);           % number of inputs nu switch transient    case 1        nu1 = nu+1;             % +1 accounts for the transient parameters    case 0        nu1 = nu;               % no transient parameters are estimatedend % switch% half the frequency width in DFT samples of the polynomial methodnn = ceil((ny + (R+1)*nu1 - 1 + mm)/2);% number of degrees of freedom in the residualsqq = 2*nn+1 - (R+1)*nu1;% equivalent number of independent experimentsM = qq + 1;% covariance matricesCY = struct('n', zeros(ny, ny, Fselect), 'm', zeros(ny, ny, Fselect), 'm_nt', zeros(ny, ny, Fselect));% sample mean outputY = struct('m', zeros(ny, Fselect), 'm_nt', zeros(ny, Fselect));% plant and noise transient contribution at outputTY = zeros(ny, Fselect);% frequency response matrixG = zeros(ny, nu, Fselect);% covariance matrix FRMCvecG = zeros(ny*nu, ny*nu, Fselect);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculation of the regressor matrix Kn %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% regressor matrixKn = zeros((R+1)*nu1, 2*nn+1);% intermediate variablePower_r = ones(R+1, 2*nn+1);% loop over all frequenciesfk = 0;                                     % frequency index of the output parametersfor kk = 1:Fstep:F        fk = fk + 1;        % range of DFT frequencies around kk    if kk <= nn        r_index = [-kk+1:1:2*nn-kk+1];    end % if    if (kk >= nn+1) & (kk <= F-nn)        r_index = [-nn:1:nn];    end % if    if kk >= F-nn+1        r_index = [-2*nn+F-kk:1:F-kk];    end % if        % intermediate variable: powers of r    r_power = data.freq(kk+r_index)-data.freq(kk);    for ii = 2:R+1        Power_r(ii,:) = (r_power).^(ii-1);    end % ii        % regressor matrix    if nu > 0        Ukr = data.U(:, kk+r_index);        for jj = 1:2*nn+1            Kn(1:(R+1)*nu, jj) = kron(Power_r(:, jj), Ukr(:, jj));        end % jj     end % if nu > 0    if transient        Kn((R+1)*nu+1:end, :) = Power_r;    end % if transient        % normalise the rows of Kn for improving the numerical stability of    % the calculations    Scale = sum(abs(Kn.^2), 2).^0.5;                            % 2-norm rows Kn    FindZeros = find(Scale == 0);                               % if the input is exactly zero in the band kk-n:kk+n    Scale(FindZeros) = 1;                                       % then the scaling is set equal to one    Kn = Kn./repmat(Scale, [1, 2*nn+1]);    % numerical stable LS estimate output (= "sample mean")    [Un, Sn, Vn] = svd(Kn', 0);    Yn = (data.Y(:, kk+r_index)*Un)*Un';                        % Yn = data.Y(:, kk+r_index)*Qn; with Qn = Un * Un'    Index_kk = find(r_index == 0);    Y.m(:, fk) = Yn(:, Index_kk);        % numerical stable LS estimate of the noise covariance matrix (= "sample covariance matrix")                          En = data.Y(:, kk+r_index) - Yn;                            % LS residuals En = data.Y(:, kk+r_index)*Pn; with Pn = I2n+1 - Qn    CY.n(:,:,fk) = (En*En')/qq;        % sample covariance of the sample mean    Qnkk = Un(Index_kk, :) * Un(Index_kk, :)';    CY.m(:,:,fk) = real(Qnkk) * CY.n(:,:,fk);                   % CYm(:,:,fk) = real(Qn(Index_kk, Index_kk)) * CYn(:,:,fk);        if transient                % LS estimate transient contribution at output        ss = diag(Sn);        IndexZeros = find(ss == 0);        ss(IndexZeros) = inf;        ss = diag(1./ss);        Theta = Yn * (Un * ss * Vn');        IndexTrans = nu*(R+1)+1;                             	% position transient parameters        TY(:, fk) = Theta(:, IndexTrans) / Scale(IndexTrans);	% denormalisation parameters         Y.m_nt = Y.m - TY;                                      % output without transient                   % covariance matrix Ym-TY        % qkk = Un * Un(Index_kk, :)';        % bmm = Un * (diag(ss) .* Vn(IndexTrans, :)') / Scale(IndexTrans);  % denormalisation parameters         % difference qkk - bmm        qkk_bmm = Un * (Un(Index_kk, :)' - (diag(ss) .* Vn(IndexTrans, :)')/Scale(IndexTrans));        CY.m_nt(:,:,fk) = norm(qkk_bmm, 2)^2 * CY.n(:,:,fk);            end % if transient        if nu > 0                % estimate FRM        IndexFRM = [1:nu];        G(:, :, fk) = Theta(:, IndexFRM) ./ repmat(Scale(IndexFRM).', [ny, 1]);        % covariance matrix vec(G)        dimVn = size(Vn, 2);                VV =  (Vn(IndexFRM, :) ./ repmat(Scale(IndexFRM), [1, dimVn])) * ss;    % intermediate variable         CvecG(:, :, fk) = kron(conj(VV * VV'), CY.n(:,:,fk));        end % nu > 0        %if kk==347;'rik 290';keyboard;end    end % kk iteration over all frequencies