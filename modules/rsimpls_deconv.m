function decon_res = rsimpls_decon(Y, X)
% RSIMPLS is a 'Robust method for Partial Least Squares Regression based on the
% SIMPLS algorithm'. Maltab code can download from https://github.com/mwgeurts/libra.git
% Optional input arguments: 
%		Y: mixture sample data
%		X: pure sample data
%	Out: deconvolution results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pvalue = calc_pval(X, L)
	uniq_label = unique(L);
	k = length(uniq_label);
	n_j = histc(L, uniq_label);
	
	sn_pos = {};
	for i = 1 : k
		sn_pos{i} = find(L == uniq_label(i));
	end

	s_j = []; x_j = [];
	for i = 1 : k
		s_j = [s_j std(X(sn_pos{i}))];
		x_j = [x_j mean(X(sn_pos{i}))];
	end

	m_j = n_j - 1;
	D = (sum(m_j * (s_j / x_j))) / sum(m_j);
	D_AD = (sum(m_j * (s_j / x_j - D)^2 )) / ( D^2 * (0.5 + D^2) );
	pvalue = chi2cdf(k - 1, D_AD, 'upper');
end

X = readtable(X, 'Delimiter', '\t', 'ReadRowNames', 1);
Y = readtable(Y, 'Delimiter', '\t', 'ReadRowNames', 1);

comNames = intersect(X.Properties.RowNames, Y.Properties.RowNames);
X = cell2mat(table2cell(X(comNames, :)));
Y = cell2mat(table2cell(Y(comNames, :)));

X_norm = sqrt(1 + X);
Y_norm =sqrt(1 + Y);

X_norm = bsxfun(@minus, X, mean(X)) ./ std(reshape(X, [], 1));
Y_norm = bsxfun(@rdivide, bsxfun(@minus, Y, mean(Y)), std(Y));

[m, n] = size(X);
results = rsimpls(X_norm, Y_norm, 'k', n, 'plots',0, 'alpha', 0.99);
results.slope = real(results.slope);
results.slope(find(results.slope < 0)) = 0;
coeffs  = (bsxfun(@rdivide, results.slope, sum(results.slope, 1)))';
Y_res   = real(results.res);
rmse    = mean(Y_res.^2);
rsquare = 1 - sum(Y_res.^2) ./ sum(bsxfun(@minus, Y_norm, mean(Y_norm)).^2);
rsq_adj = 1 - ((1 - rsquare) * (m - 1) / (m - n - 1));

pvalues = [];
L = strcat(repmat('a', 1, size(Y, 1)), repmat('b', 1, size(Y, 1))); 
results.fitted = real(results.fitted);

for i = 1 : size(Y, 2)
	vec = [ Y_norm(:, i); results.fitted(:, i) ]';
	pvalues = [ pvalues, calc_pval(vec, L) ];
end	
decon_res = [ coeffs, rsq_adj', rmse', pvalues' ];
end

