function [z] = lognormpdf(x, mu, sigma2)
    z = -0.5*log(2*pi*sigma2) -(x-mu).^2./(2*sigma2);
end