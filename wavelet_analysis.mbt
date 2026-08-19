///|
/// Continuous wavelet families supported by the portable analyzer.
pub enum WaveletFamily {
  RickerWavelet
  MorletWavelet
  GaussianDerivative
}

///|
/// Return wavelet family tags.
pub fn wavelet_ricker() -> WaveletFamily {
  RickerWavelet
}

///|
pub fn wavelet_morlet() -> WaveletFamily {
  MorletWavelet
}

///|
pub fn wavelet_gaussian_derivative() -> WaveletFamily {
  GaussianDerivative
}

///|
/// One coefficient in a time-scale representation.
pub struct WaveletCoefficient {
  sample_index : Int
  scale : Double
  value : Double
  magnitude : Double
}

///|
/// A compact scalogram with row-major coefficients.
pub struct Scalogram {
  scales : Array[Double]
  sample_count : Int
  coefficients : Array[Double]
  family : WaveletFamily
}

///|
/// Evaluate a selected wavelet at a normalized time coordinate.
fn wavelet_value(family : WaveletFamily, time : Double) -> Double {
  match family {
    RickerWavelet => {
      let squared = time * time
      (1.0 - 2.0 * squared) * @math.exp(-squared)
    }
    MorletWavelet => @math.cos(5.0 * time) * @math.exp(-0.5 * time * time)
    GaussianDerivative => -time * @math.exp(-0.5 * time * time)
  }
}

///|
/// Compute a real-valued continuous wavelet transform.
pub fn continuous_wavelet_transform(
  data : Array[Double],
  scales : Array[Double],
  family : WaveletFamily,
) -> Scalogram {
  let coefficients = Array::make(data.length() * scales.length(), 0.0)
  for scale_index = 0
      scale_index < scales.length()
      scale_index = scale_index + 1 {
    let scale = if scales[scale_index] <= 0.0 {
      1.0
    } else {
      scales[scale_index]
    }
    let radius = (3.0 * scale).to_int()
    for sample = 0; sample < data.length(); sample = sample + 1 {
      let mut total = 0.0
      let mut weight = 0.0
      for offset = -radius; offset <= radius; offset = offset + 1 {
        let source = sample + offset
        if source >= 0 && source < data.length() {
          let argument = offset.to_double() / scale
          let value = wavelet_value(family, argument) / scale.sqrt()
          total += data[source] * value
          weight += value * value
        }
      }
      coefficients[scale_index * data.length() + sample] = if weight == 0.0 {
        0.0
      } else {
        total / weight.sqrt()
      }
    }
  }
  { scales, sample_count: data.length(), coefficients, family }
}

///|
/// Read one coefficient from a scalogram, returning zero when out of bounds.
pub fn scalogram_at(
  scalogram : Scalogram,
  scale_index : Int,
  sample_index : Int,
) -> Double {
  if scale_index < 0 ||
    scale_index >= scalogram.scales.length() ||
    sample_index < 0 ||
    sample_index >= scalogram.sample_count {
    0.0
  } else {
    scalogram.coefficients[scale_index * scalogram.sample_count + sample_index]
  }
}

///|
/// Return energy by wavelet scale.
pub fn wavelet_scale_energy(scalogram : Scalogram) -> Array[Double] {
  let result = Array::make(scalogram.scales.length(), 0.0)
  for scale_index = 0
      scale_index < scalogram.scales.length()
      scale_index = scale_index + 1 {
    let mut energy = 0.0
    for sample = 0; sample < scalogram.sample_count; sample = sample + 1 {
      let value = scalogram_at(scalogram, scale_index, sample)
      energy += value * value
    }
    result[scale_index] = energy
  }
  result
}

///|
/// Return the scale with the greatest transform energy.
pub fn dominant_wavelet_scale(scalogram : Scalogram) -> Double {
  if scalogram.scales.is_empty() {
    return 0.0
  }
  let energy = wavelet_scale_energy(scalogram)
  let mut best = 0
  for index = 1; index < energy.length(); index = index + 1 {
    if energy[index] > energy[best] {
      best = index
    }
  }
  scalogram.scales[best]
}

///|
/// Find the largest wavelet coefficient and its location.
pub fn strongest_wavelet_coefficient(
  scalogram : Scalogram,
) -> WaveletCoefficient {
  let mut best_scale = 0
  let mut best_sample = 0
  let mut best_value = 0.0
  for scale_index = 0
      scale_index < scalogram.scales.length()
      scale_index = scale_index + 1 {
    for sample = 0; sample < scalogram.sample_count; sample = sample + 1 {
      let value = scalogram_at(scalogram, scale_index, sample)
      if value.abs() > best_value.abs() {
        best_value = value
        best_scale = scale_index
        best_sample = sample
      }
    }
  }
  {
    sample_index: best_sample,
    scale: if scalogram.scales.is_empty() {
      0.0
    } else {
      scalogram.scales[best_scale]
    },
    value: best_value,
    magnitude: best_value.abs(),
  }
}

///|
/// Soft-threshold transform coefficients and reconstruct a denoised trace.
pub fn wavelet_denoise(
  data : Array[Double],
  scales : Array[Double],
  family : WaveletFamily,
  threshold : Double,
) -> Array[Double] {
  let scalogram = continuous_wavelet_transform(data, scales, family)
  let result = Array::make(data.length(), 0.0)
  let limit = if threshold < 0.0 { 0.0 } else { threshold }
  if scales.is_empty() {
    return result
  }
  for sample = 0; sample < data.length(); sample = sample + 1 {
    let mut total = 0.0
    for scale_index = 0
        scale_index < scales.length()
        scale_index = scale_index + 1 {
      let coefficient = scalogram_at(scalogram, scale_index, sample)
      let magnitude = coefficient.abs()
      let shrink = if magnitude <= limit {
        0.0
      } else {
        coefficient.signum() * (magnitude - limit)
      }
      total += shrink
    }
    result[sample] = total / scales.length().to_double()
  }
  result
}

///|
/// Build logarithmically spaced scales for multi-resolution analysis.
pub fn logarithmic_scales(
  minimum : Double,
  maximum : Double,
  count : Int,
) -> Array[Double] {
  let result = Array::new()
  if minimum <= 0.0 || maximum < minimum || count <= 0 {
    return result
  }
  if count == 1 {
    result.push(minimum)
    return result
  }
  let ratio = @math.ln(maximum / minimum) / (count - 1).to_double()
  for index = 0; index < count; index = index + 1 {
    result.push(minimum * @math.exp(ratio * index.to_double()))
  }
  result
}
