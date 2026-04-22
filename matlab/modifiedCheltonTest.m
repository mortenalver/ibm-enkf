

N = size(values,1);
X = values(:,1);
Y = values(:,2);
% 2. Calculate autocorrelation (up to a sufficient lag, e.g., N/2 or N/4)
% It is common to use autocorrelation up to where it becomes negligible
[acfX, lagsX] = autocorr(X, 'NumLags', N-1);
[acfY, lagsY] = autocorr(Y, 'NumLags', N-1);

% 3. Calculate Effective Sample Size (N*)
% Using the formula: N / sum(rhoX * rhoY)
% NOTE: For many, you might sum up to a specific lag limit (e.g., Pyper & Peterman)
prodAcf = acfX .* acfY;
effectiveN = N / sum(prodAcf);

% 4. Compute Sample Correlation (r)
[r, p_original] = corr(X, Y);

% 5. Compute Corrected Significance (p-value)
% t = r * sqrt((N* - 2) / (1 - r^2))
t_stat = r * sqrt((effectiveN - 2) / (1 - r^2));
p_corrected = 2 * (1 - tcdf(abs(t_stat), effectiveN - 2));

% Output
fprintf('Original Sample Size: %d\n', N);
fprintf('Effective Sample Size (N*): %.2f\n', effectiveN);
fprintf('Correlation Coefficient: %.3f\n', r);
fprintf('Corrected P-value: %.5f\n', p_corrected);