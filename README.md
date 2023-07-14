# sgolay_robust
A Savitzky-Golay filter with robust options

`sgolay_robust` is a **MATLAB** script that permits to smooth data thanks to local polynomial regressions based on the **iteratively-reweighted least square** (IRLS) procedure. The IRLS method minimizes the influence of outliers using *Bisquare weighting*.

The filter can be used for various applications, such as anomaly detection, baseline correction, peak detection...

The general syntax is:

`[output_1, output_2] = sgolay_robust(input_data,frame_half_length,order,evaluation_interval,options)`

where:

* `input_data` is a vector of real numbers (*i.e.* the data to smooth.)
* `frame_half_length` is a scalar positive integer that specifies the half-length of the filter frame. If the value equals 10, the filter will operate on each value of `input_data` using 21 values (the 10 previous values, the 10 subsequent values and the value itself).
* `order` is a scalar positive or null integer that specifies the polynomial order for the filter.
* `evaluation_interval` is a scalar positive integer that specifies the interval for which the filter operates. When set to 1, the filter smoothes the data on every value. When set to 10, the filter smoothes the data every 10 points, and evaluates the intermediate values by interpolation.
* `options`: optional name-values arguments
  * `interpolation_method`: a scalar string that specifies the interpolation method, in case `evaluation_interval` is larger than 1. The possible interpolation methods are `"linear"` (default), `"cubic"`, `"spline"` and `"makima"`.
  * `display_warning`: a scalar string that specifies whether or not the console should display warning messages.
  * `TolFun`: a scalar positive double that specifies the stopping threshold in the change of the objective function value during the IRLS procedure. The default value is `1e-3`.
  * `TolX`: a scalar positive double that specifies the stopping threshold in the change of the polynomial coefficients during the IRLS procedure. The default value is `1e-3`.
  * `MaxIter`: a scalar positive integer that specifies the maximum allowed iterations for the IRLS procedure. The default value is `200`.
 
The two returned outputs are:
* `output_1`: the smoothed data
* `output_2`: the local robust standard deviation of the smoothed data

Examples of function calls:
* `[y0, s] = sgolay_robust(y,10,2,1)`
* `[y0, s] = sgolay_robust(y,10,2,10,"interpolation_method","makima")`
* `[y0, s] = sgolay_robust(y,10,2,10,"interpolation_method","makima","MaxIter",50)`

## Warning

The parameters of the filter should be appropriately set to avoid overfitting. Using a large polynomial order value is therefore not a good idea, even if the script performs a certain regularization of the data during the IRLS procedure (matrix inversion according to Moore-Penrose).

## Examples

Examples are provided [here](./EXAMPLES.md).

## Disclaimer

Author: Antoine BUREL

The script is licensed under the CC0 1.0 Universal terms. See the [LICENSE](./LICENSE) file.

The author accepts no responsibility for the use of the script.
