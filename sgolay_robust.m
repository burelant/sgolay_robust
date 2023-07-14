function [y1,y2,y3] = sgolay_robust(input_data,frame_half_length,order,evaluation_interval,opts)

% SGOLAY_ROBUST  Smooth data using Savitzky-Golay filter with
% iteratively-reweighted least-square procedure implementation to reduce
% the influence of input data outliers.
%
%   [y,s,p] = sgolay_robust(input_data,frame_half_length,order,evaluation_interval,opts) returns a
%   vector y containing input_data smoothed data, a vector s containing the
%   local residual standard errors, and a matrix p containing the
%   polynomial coefficients for each point. The parameters of the filter
%   are:
%   - frame_half_length: the length of the unilateral frame of the filter.
%   - order: the polynomial order.
%   - evaluation_interval: the evaluation_interval size according to which the data are smoothed. If
%   evaluation_interval is larger than unity, intermediate points are estimated by
%   interpolation.
%   - opts: a structure or a list of paired name-value arguments:
%       * display_warnings: ("on" or "off") permits to activate or
%       deactivate warnings in the console. Default value: "on".
%       * TolFun: (1e-3 by default) the tolerance value for the local
%       weighted standard error which stops the IRLS procedure.
%       * TolX: (1e-3 by default) the tolerance value for the local
%       polynomial coefficients which stops the IRLS procedure.
%       * MaxIter: (200 by default) the maximum number of iterations after
%       which the IRLS procedure stops.
%       * interpolation_method: ("linear","cubic","spline") Determines the way
%       according to which interpolations are performed, if evaluation_interval > 1.
%       Default value: "linear".

arguments % Argument validation
    input_data double {mustBeNumeric,mustBeReal,mustBeVector} 
    frame_half_length (1,1) double {mustBeNumeric,mustBeReal,mustBePositive}
    order (1,1) double {mustBeNumeric,mustBeReal, ...
        mustBeGreaterThanOrEqual(order,0)}
    evaluation_interval (1,1) double {mustBeNumeric,mustBeReal,mustBePositive} 
    opts.interpolation_method string {mustBeTextScalar,mustBeMember( ...
        opts.interpolation_method,["linear","cubic","spline","makima"])} = "linear"
    % opts.RobustWeightFunction string {mustBeTextScalar,mustBeMember(opts.RobustWeightFunction,["bisquare"])} = "bisquare"
    opts.display_warnings string {mustBeTextScalar, ...
        mustBeMember(opts.display_warnings,["on","off"])} = "on"
    opts.TolFun (1,1) double ...
        {mustBeReal,mustBePositive,mustBeNumeric} = 1e-3 ;
    opts.TolX (1,1) double ...
        {mustBeReal,mustBePositive,mustBeNumeric} = 1e-3 ;
    opts.MaxIter (1,1) double ...
        {mustBeReal,mustBePositive,mustBeNumeric} = 200 ;
end

% Disable warnings if requested
if opts.display_warnings == "off"
    warning("off","all") ;
end
% When appropriate, w is rounded to nearest integer
if mod(frame_half_length,1) ~= 0
    frame_half_length = round(frame_half_length) ;
    warning("Frame half-length was rounded to nearest integer.") ;
end
% When appropriate, o is rounded to nearest integer
if mod(order,1) ~= 0
    order = round(order) ;
    warning("Filter order was rounded to nearest integer.") ;
end
% When appropriate, evaluation_interval is rounded to nearest integer
if mod(evaluation_interval,1) ~= 0
    evaluation_interval = min(1,round(evaluation_interval)) ;
    warning("evaluation_interval size was rounded to nearest integer.") ;
end

% Ckeck for NaNs
test = sum(isnan(input_data)) ;
if test > 0
    error("The data supplied contain NaNs. Operation aborted.") ;
    return
end

test = isrow(input_data) ;
if test == true
    input_data_type = "row" ;
    input_data = input_data' ;
else
    input_data_type = "column" ;
end

% Memory preallocation
y1 = zeros(size(input_data)) ;
y2 = zeros(size(input_data)) ;
y3 = zeros(numel(input_data),order+1) ;

N = numel(input_data) ;

for i = 1:evaluation_interval:N

    ind1 = max(1,i-frame_half_length) ;
    ind2 = min(N,i+frame_half_length) ;
    y_spl = input_data(ind1:ind2) ;
    x_spl = [ind1:ind2]' ;

    weights = ones(size(y_spl)) ;

    p = inf(order+1,1) ; s = inf ; ds = 1 ; dp = 1 ;

    X = repmat(x_spl,1,order+1).^[order:-1:0] ;
    c = 4.685 ; nb_iter = 0 ;

    while (abs(ds) >= opts.TolFun) || (abs(dp) >= opts.TolX)

        if nb_iter > opts.MaxIter
            break
        end

        nb_iter = nb_iter+1 ;

        old_p = p ;
        p = pinv(sqrt(weights).*X)*(sqrt(weights).*y_spl) ;
    
        y_spl_calc = X*p ;
        r = y_spl_calc-y_spl ;
        tau = median(abs(r-median(r)))/0.6745 ;
        if tau == 0
            s = std(r)*sqrt(1+1/numel(y_spl)+(i-mean(x_spl)).^2/sum((x_spl-mean(x_spl)).^2)) ;
            break
        end
        z = r/tau ;
        weights = (abs(z)<c).*(1-(z/c).^2).^2 ;
        weights = weights.*(weights>=0) ; 
        old_s = s ;
        s = std(r,weights)*sqrt(1+1/numel(y_spl)+(i-sum(weights.*x_spl)./sum(weights)).^2/sum((x_spl-sum(weights.*x_spl)./sum(weights)).^2)) ;
        dp = max((old_p-p)./p) ;
        ds = (old_s-s)/s ;

    end

    y3(i,:) = p' ;
    y2(i,1) = s ;
    y1(i,1) = polyval(p,i) ;
 
end

if input_data_type == "row"
    y1 = y1' ;
    y2 = y2' ;
    x3 = y3' ;
    input_data = input_data' ;
    x_interp = 1:evaluation_interval:N ;
    x_query = 1:N ;
else
    x_interp = [1:evaluation_interval:N]' ;
    x_query = [1:N]' ;
end

if evaluation_interval == 1
    % No need for interpolation.
    return
end

y1 = interp1(x_interp,y1(1:evaluation_interval:N),x_query,opts.interpolation_method,"extrap") ;
y2 = interp1(x_interp,y2(1:evaluation_interval:N),x_query,opts.interpolation_method,"extrap") ;
y3 = interp1(x_interp,y3(1:evaluation_interval:N,:),x_query,opts.interpolation_method,"extrap") ;

end