function Value = goodnessOfFit_edit(X,Xref,Measure)
%GOODNESSOFFIT Computes a goodness of fit measure for data.
%
%  FIT = GOODNESSOFFIT(X, XREF, MEASURE)
%    compares the data matrix X with a reference value XREF using MEASURE
%    as a comparison yardstick to produce the number FIT. FIT is a
%    quantitative representation of the closeness of X to XREF.
%
%  Input arguments:
%           X: Test data matrix with Ns samples and N channels; Ns-by-N
%              matrix. X must not contain NaNs or Infs.
%        XREF: Reference data of same size as X. Must not contain NaNs or
%              Infs.
%     MEASURE: One of 'MSE', 'NRMSE', NMSE'. The formulas used for these
%              measures are as follows:
%              'MSE': FIT = norm(X-XREF)^2/Ns;
%            'NRMSE': FIT(i) = 1 - (norm(XREF(:,i)-X(:,i)))/(norm(XREF(:,i)-mean(XREF(:,i))))
%             'NMSE': FIT(i) = 1 - (norm(XREF(:,i)-X(:,i)))^2/(norm(XREF(:,i)-mean(XREF(:,i))))^2
%
%    The output argument FIT is a double scalar if MEASURE = 'MSE', and a
%    row vector with N entries for other values of MEASURE, where N =
%    size(X,2). NMSE and NRMSE vary between -Inf (bad fit) to 1 (perfect
%    fit); a zero value for NMSE or NRMSE denotes that the data X is no
%    better than a straight line at matching XREF. 'NRMSE' is the measure
%    employed by the COMPARE function.
%
%  Comparing multiple data sets:
%    Either X or XREF or both can be a cell array of double matrices if
%    you want to compare one or more test data sets to one or more
%    reference values. In that case, FIT is a double array of size equal to
%    [Nx Cell_Array_Size], where Nx is 1 for MEASURE = 'MSE' and N for
%    other values of MEASURE. If both XREF are cell arrays, their sizes
%    must match so that FIT(:,i1,i2,...) corresponds to X{i1,i2,...} and
%    XREF{i1, i2,...}.
%
% See also COMPARE, PE, RESID, FPE, AIC.

%   Author(s): Lennart Ljung, Rajiv Singh
%   Copyright 1986-2016 The MathWorks, Inc.

if ~iscell(Xref), Xref = {Xref}; end
if ~iscell(X), X = {X}; end

if ~isequal(size(X),size(Xref)) && ~isscalar(X) && ~isscalar(Xref)
   ctrlMsgUtils.error('Ident:analysis:goodnessoffit1')
else
   for ct = 1:numel(X)
      if ~isequal(size(X{min(ct,1)}),size(Xref{min(ct,1)}))
         ctrlMsgUtils.error('Ident:analysis:goodnessoffit2')
      end
   end
end

Ne = length(Xref);
if strcmpi(Measure,'mse')
   Value = zeros(1,Ne);
else
   Value = zeros(size(Xref{1},2),Ne); % N-by-Ne
end

for ct = 1:Ne
   x1 = X{min(ct,end)}; x2 = Xref{min(ct,end)};
   N = size(x1,2);
   switch lower(Measure)
      case 'mse'
         e = x1-x2;
         Value(ct) = trace(e'*e)/(size(x1,1));
      case {'nrmse','nmse'}
         for k = 1:N
            n1 = norm(x2(:,k)-x1(:,k),'fro');
%             n2 = norm(x2(:,k)-mean(x2(:,k)));
            n2 = norm(x2(:,k),'fro');
            if n2==0 && n1==0
               Rat = 0;
            else
               Rat = n1/n2;
            end
            if strcmpi(Measure,'nmse')
               Rat = Rat^2;
            end
            Value(k,ct) = 1 - Rat;
         end
   end
end
