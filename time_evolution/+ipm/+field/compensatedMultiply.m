function [high,low] = compensatedMultiply(A,U,Ulow)
%IPM.FIELD.COMPENSATEDMULTIPLY Sparse stencil times a double-double field.
%   [HIGH,LOW] = compensatedMultiply(A,U,ULOW) represents A*(U+ULOW)
%   without rounding away the trailing component of the field. A is sparse
%   real double, and U is dense real double. ULOW is optional and denotes
%   the normalized trailing part of U, normally obtained with TwoSum.
%   No coefficients or diagonal entries are replaced. Products use Dekker
%   splitting; the stencil sum uses a double-double accumulator. Sparse
%   diagonals are evaluated in vectorized batches rather than node loops.
%
%   For a right stencil U*A.', transpose both parts of U, call this helper,
%   then transpose both returned parts. Keep HIGH and LOW separate when a
%   later derivative or residual must include the trailing information.

if ~isa(A,'double') || ~isreal(A) || ~issparse(A) || ...
        ~isa(U,'double') || ~isreal(U) || issparse(U) || ...
        ~ismatrix(U) || size(A,2)~=size(U,1)
    error('ipm:CompensatedMultiplyDimensions', ...
        'Expected sparse real double A and matching dense real double U.');
end
hasLow=nargin>=3 && ~isempty(Ulow);
if hasLow && (~isa(Ulow,'double') || ~isreal(Ulow) || ...
        issparse(Ulow) || ~isequal(size(Ulow),size(U)))
    error('ipm:CompensatedMultiplyDimensions', ...
        'Ulow must be a dense real double array with the same size as U.');
end
if any(~isfinite(U(:))) || (hasLow && any(~isfinite(Ulow(:))))
    error('ipm:CompensatedMultiplyFinite', ...
        'The complete input field must be finite.');
end
[rows,columns,values]=find(A);
if any(~isfinite(values))
    error('ipm:CompensatedMultiplyFinite', ...
        'The original stencil coefficients must be finite.');
end
splitter=134217729; % 2^27+1, for IEEE binary64.
safeMaximum=realmax('double')/splitter;
if any(abs(U(:))>safeMaximum) || any(abs(values)>safeMaximum)
    error('ipm:CompensatedMultiplyRange', ...
        'An operand exceeds the safe range of the double-double splitter.');
end
[uHigh,uSplitLow]=split(U,splitter);
high=zeros(size(A,1),size(U,2));
low=zeros(size(high));
offsets=columns-rows;
for offset=reshape(unique(offsets),1,[])
    selected=offsets==offset;
    destination=rows(selected);
    source=columns(selected);
    coefficient=values(selected);
    [coefficientHigh,coefficientLow]=split(coefficient,splitter);
    value=U(source,:);
    valueHigh=uHigh(source,:);
    valueLow=uSplitLow(source,:);
    product=coefficient.*value;
    productError=((coefficientHigh.*valueHigh-product)+ ...
        coefficientHigh.*valueLow+coefficientLow.*valueHigh)+ ...
        coefficientLow.*valueLow;
    if hasLow
        % Ulow is the normalized low component: the rounding error of
        % this product is already of second order relative to A*U.
        productError=productError+coefficient.*Ulow(source,:);
        [product,productError]=two_sum(product,productError);
    end
    [high(destination,:),low(destination,:)]=add_pair( ...
        high(destination,:),low(destination,:),product,productError);
end
if any(~isfinite(high(:))) || any(~isfinite(low(:)))
    error('ipm:CompensatedMultiplyFinite', ...
        'The complete stencil product is not finite.');
end
end

function [high,low]=split(value,splitter)
temporary=splitter*value;
high=temporary-(temporary-value);
low=value-high;
end

function [high,low]=add_pair(aHigh,aLow,bHigh,bLow)
% General TwoSum remains valid during cancellation, including a zero
% leading component; no magnitude-order assumption is needed.
[leading,leadingError]=two_sum(aHigh,bHigh);
[trailing,trailingError]=two_sum(aLow,bLow);
[middle,middleError]=two_sum(leadingError,trailing);
[head,headError]=two_sum(leading,middle);
tail=(headError+middleError)+trailingError;
[high,low]=two_sum(head,tail);
end

function [sum,error]=two_sum(a,b)
sum=a+b;
virtualB=sum-a;
error=(a-(sum-virtualB))+(b-virtualB);
end
