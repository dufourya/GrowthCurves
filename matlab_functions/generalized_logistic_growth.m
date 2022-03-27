function OD = generalized_logistic_growth(t, A, K, C, Q, B, v)

OD = A+(K-A)./((C+exp(-t/B)./Q).^(1/v));

end