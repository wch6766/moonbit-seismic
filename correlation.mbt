///| Result of a cross-correlation operation.
pub struct CorrelationResult {
  lags : Array[Int]
  coefficients : Array[Double]
  max_coefficient : Double
  max_lag : Int
  max_time_offset : Double
}

///| Calculate standard cross-correlation between two signals.
/// Returns an array of correlation values for each lag from -N to +N, where N is max_lag.
pub fn cross_correlate(
  signal : Array[Double],
  template : Array[Double],
  max_lag : Int
) -> CorrelationResult {
  let sig_len = signal.length()
  let temp_len = template.length()
  if sig_len == 0 || temp_len == 0 || max_lag < 0 {
    return {
      lags: [],
      coefficients: [],
      max_coefficient: 0.0,
      max_lag: 0,
      max_time_offset: 0.0,
    }
  }
  
  let num_lags = 2 * max_lag + 1
  let lags = Array::new(capacity=num_lags)
  let coefficients = Array::new(capacity=num_lags)
  
  let mut max_coeff = -2.0
  let mut best_lag = 0
  
  // Calculate template stats once
  let mut temp_sum = 0.0
  let mut temp_sq_sum = 0.0
  for v in template {
    temp_sum = temp_sum + v
    temp_sq_sum = temp_sq_sum + v * v
  }
  let temp_mean = temp_sum / temp_len.to_double()
  let temp_variance = temp_sq_sum - temp_sum * temp_mean
  let temp_std = if temp_variance <= 0.0 { 0.0 } else { temp_variance.sqrt() }
  
  for lag = -max_lag; lag <= max_lag; lag = lag + 1 {
    lags.push(lag)
    
    // We compute the normalized cross-correlation (NCC) at this lag
    let mut dot_product = 0.0
    let mut sig_sum = 0.0
    let mut sig_sq_sum = 0.0
    let mut overlap_count = 0
    
    for i = 0; i < temp_len; i = i + 1 {
      let sig_idx = i + lag
      if sig_idx >= 0 && sig_idx < sig_len {
        let x = signal[sig_idx]
        let y = template[i]
        dot_product = dot_product + x * y
        sig_sum = sig_sum + x
        sig_sq_sum = sig_sq_sum + x * x
        overlap_count = overlap_count + 1
      }
    }
    
    let coeff = if overlap_count == 0 || temp_std == 0.0 {
      0.0
    } else {
      let sig_mean = sig_sum / overlap_count.to_double()
      let sig_variance = sig_sq_sum - sig_sum * sig_mean
      let sig_std = if sig_variance <= 0.0 { 0.0 } else { sig_variance.sqrt() }
      
      if sig_std == 0.0 {
        0.0
      } else {
        // Normalized cross correlation
        let mut cov = 0.0
        for i = 0; i < temp_len; i = i + 1 {
          let sig_idx = i + lag
          if sig_idx >= 0 && sig_idx < sig_len {
            cov = cov + (signal[sig_idx] - sig_mean) * (template[i] - temp_mean)
          }
        }
        cov / (sig_std * temp_std)
      }
    }
    
    coefficients.push(coeff)
    if coeff > max_coeff {
      max_coeff = coeff
      best_lag = lag
    }
  }
  
  {
    lags,
    coefficients,
    max_coefficient: max_coeff,
    max_lag: best_lag,
    max_time_offset: 0.0, // Calculated by caller who knows sampling rate
  }
}

///| Slide a template trace across a main trace to compute sliding-window cross correlation.
pub fn correlate_traces(trace : Trace, template : Trace) -> CorrelationResult {
  let max_lag = trace.npts() - 1
  let res = cross_correlate(trace.data, template.data, max_lag)
  
  // Calculate exact time offset for the best lag
  let dt = trace.delta()
  let max_time_offset = res.max_lag.to_double() * dt + (trace.start_time - template.start_time)
  
  {
    lags: res.lags,
    coefficients: res.coefficients,
    max_coefficient: res.max_coefficient,
    max_lag: res.max_lag,
    max_time_offset,
  }
}
