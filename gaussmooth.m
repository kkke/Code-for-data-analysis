%% temporal filting the trace with a gaussian filter
function fdata = gaussmooth(data,size,sigma)
if nargin <2
    sigma = 2;
    size = 5;    % length of gaussFilter vector
end
x = linspace(-size / 2, size / 2, size);
gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
gaussFilter = gaussFilter / sum (gaussFilter); % normalize
fdata       = conv(data,gaussFilter,'same');

