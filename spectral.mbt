///|
/// A small complex number type used by the portable spectral routines.
pub struct ComplexNumber {
  real : Double
  imaginary : Double
}

///|
/// Construct a complex number from Cartesian components.
pub fn ComplexNumber::new(real : Double, imaginary : Double) -> ComplexNumber {
  { real, imaginary }
}

///|
/// Add two complex values.
pub fn ComplexNumber::add(
  self : ComplexNumber,
  other : ComplexNumber,
) -> ComplexNumber {
  { real: self.real + other.real, imaginary: self.imaginary + other.imaginary }
}

///|
/// Multiply two complex values.
pub fn ComplexNumber::mul(
  self : ComplexNumber,
  other : ComplexNumber,
) -> ComplexNumber {
  {
    real: self.real * other.real - self.imaginary * other.imaginary,
    imaginary: self.real * other.imaginary + self.imaginary * other.real,
  }
}

///|
/// Return the complex conjugate.
pub fn ComplexNumber::conjugate(self : ComplexNumber) -> ComplexNumber {
  { real: self.real, imaginary: -self.imaginary }
}

///|
/// Return the squared magnitude.
pub fn ComplexNumber::magnitude_squared(self : ComplexNumber) -> Double {
  self.real * self.real + self.imaginary * self.imaginary
}

///|
/// Return the magnitude.
pub fn ComplexNumber::magnitude(self : ComplexNumber) -> Double {
  self.magnitude_squared().sqrt()
}

///|
/// Return the phase angle in radians.
pub fn ComplexNumber::phase(self : ComplexNumber) -> Double {
  @math.atan2(self.imaginary, self.real)
}

///|
/// One-sided spectrum with frequency and power arrays.
pub struct Spectrum {
  frequencies : Array[Double]
  amplitudes : Array[Double]
  powers : Array[Double]
  phases : Array[Double]
  sample_rate : Double
  input_length : Int
}

///|
/// Return the positive-frequency DFT bins for a real-valued signal.
///
/// This implementation is intentionally backend-portable and has no native
/// FFI dependency. For long traces, callers can use windows or decimation
/// before requesting a spectrum.
pub fn real_dft(data : Array[Double]) -> Array[ComplexNumber] {
  let n = data.length()
  let count = n / 2 + 1
  let result = Array::new(capacity=count)
  if n == 0 {
    return result
  }
  for frequency_bin = 0
      frequency_bin < count
      frequency_bin = frequency_bin + 1 {
    let mut real = 0.0
    let mut imaginary = 0.0
    for sample = 0; sample < n; sample = sample + 1 {
      let angle = 2.0 *
        @math.PI *
        frequency_bin.to_double() *
        sample.to_double() /
        n.to_double()
      real += data[sample] * @math.cos(angle)
      imaginary -= data[sample] * @math.sin(angle)
    }
    result.push({ real, imaginary })
  }
  result
}

///|
/// Return a full complex DFT, including negative frequencies.
pub fn complex_dft(data : Array[ComplexNumber]) -> Array[ComplexNumber] {
  let n = data.length()
  let result = Array::new(capacity=n)
  if n == 0 {
    return result
  }
  for frequency_bin = 0; frequency_bin < n; frequency_bin = frequency_bin + 1 {
    let mut sum = { real: 0.0, imaginary: 0.0 }
    for sample = 0; sample < n; sample = sample + 1 {
      let angle = 2.0 *
        @math.PI *
        frequency_bin.to_double() *
        sample.to_double() /
        n.to_double()
      let rotation = { real: @math.cos(angle), imaginary: -@math.sin(angle) }
      sum = sum.add(data[sample].mul(rotation))
    }
    result.push(sum)
  }
  result
}

///|
/// Compute an inverse DFT for a full complex spectrum.
pub fn inverse_complex_dft(data : Array[ComplexNumber]) -> Array[ComplexNumber] {
  let n = data.length()
  let result = Array::new(capacity=n)
  if n == 0 {
    return result
  }
  for sample = 0; sample < n; sample = sample + 1 {
    let mut sum = { real: 0.0, imaginary: 0.0 }
    for frequency_bin = 0; frequency_bin < n; frequency_bin = frequency_bin + 1 {
      let angle = 2.0 *
        @math.PI *
        frequency_bin.to_double() *
        sample.to_double() /
        n.to_double()
      let rotation = { real: @math.cos(angle), imaginary: @math.sin(angle) }
      sum = sum.add(data[frequency_bin].mul(rotation))
    }
    result.push({
      real: sum.real / n.to_double(),
      imaginary: sum.imaginary / n.to_double(),
    })
  }
  result
}

///|
/// Build a one-sided periodogram with a selected analysis window.
pub fn periodogram(
  data : Array[Double],
  sample_rate : Double,
  window : Array[Double],
) -> Spectrum {
  let count = data.length()
  let safe_rate = if sample_rate > 0.0 { sample_rate } else { 1.0 }
  let tapered = data.copy()
  if window.length() == tapered.length() {
    apply_window(tapered, window)
  }
  let bins = real_dft(tapered)
  let frequencies = Array::new(capacity=bins.length())
  let amplitudes = Array::new(capacity=bins.length())
  let powers = Array::new(capacity=bins.length())
  let phases = Array::new(capacity=bins.length())
  for index = 0; index < bins.length(); index = index + 1 {
    let bin = bins[index]
    let scale = if count == 0 { 0.0 } else { 2.0 / count.to_double() }
    let amplitude = bin.magnitude() * scale
    frequencies.push(index.to_double() * safe_rate / count.to_double())
    amplitudes.push(amplitude)
    powers.push(bin.magnitude_squared() / count.to_double())
    phases.push(bin.phase())
  }
  {
    frequencies,
    amplitudes,
    powers,
    phases,
    sample_rate: safe_rate,
    input_length: count,
  }
}

///|
/// Compute a Hann-windowed periodogram.
pub fn hann_periodogram(data : Array[Double], sample_rate : Double) -> Spectrum {
  periodogram(data, sample_rate, hann_window(data.length()))
}

///|
/// Return the frequency with the largest non-DC amplitude.
pub fn dominant_frequency(spectrum : Spectrum) -> Double {
  if spectrum.amplitudes.length() < 2 {
    return 0.0
  }
  let mut best = 1
  for index = 2; index < spectrum.amplitudes.length(); index = index + 1 {
    if spectrum.amplitudes[index] > spectrum.amplitudes[best] {
      best = index
    }
  }
  spectrum.frequencies[best]
}

///|
/// Integrate spectral power over a frequency interval.
pub fn band_power(
  spectrum : Spectrum,
  low_frequency : Double,
  high_frequency : Double,
) -> Double {
  let low = if low_frequency < 0.0 { 0.0 } else { low_frequency }
  let high = if high_frequency < low { low } else { high_frequency }
  let mut total = 0.0
  for index = 0; index < spectrum.powers.length(); index = index + 1 {
    let frequency = spectrum.frequencies[index]
    if frequency >= low && frequency <= high {
      total += spectrum.powers[index]
    }
  }
  total
}

///|
/// Return the power-weighted spectral centroid.
pub fn spectral_centroid(spectrum : Spectrum) -> Double {
  let mut weighted = 0.0
  let mut total = 0.0
  for index = 0; index < spectrum.powers.length(); index = index + 1 {
    let power = spectrum.powers[index]
    weighted += spectrum.frequencies[index] * power
    total += power
  }
  if total == 0.0 {
    0.0
  } else {
    weighted / total
  }
}

///|
/// Return the power-weighted spectral spread around the centroid.
pub fn spectral_spread(spectrum : Spectrum) -> Double {
  let centroid = spectral_centroid(spectrum)
  let mut weighted = 0.0
  let mut total = 0.0
  for index = 0; index < spectrum.powers.length(); index = index + 1 {
    let distance = spectrum.frequencies[index] - centroid
    weighted += distance * distance * spectrum.powers[index]
    total += spectrum.powers[index]
  }
  if total == 0.0 {
    0.0
  } else {
    (weighted / total).sqrt()
  }
}

///|
/// Return the normalized Shannon entropy of the power spectrum.
pub fn spectral_entropy(spectrum : Spectrum) -> Double {
  let total = band_power(spectrum, 0.0, spectrum.sample_rate / 2.0)
  if total == 0.0 {
    return 0.0
  }
  let mut entropy = 0.0
  for power in spectrum.powers {
    if power > 0.0 {
      let probability = power / total
      entropy -= probability * @math.ln(probability)
    }
  }
  let bins = spectrum.powers.length()
  if bins <= 1 {
    0.0
  } else {
    entropy / @math.ln(bins.to_double())
  }
}

///|
/// Return the frequency below which `fraction` of the total power lies.
pub fn spectral_rolloff(spectrum : Spectrum, fraction : Double) -> Double {
  let total = band_power(spectrum, 0.0, spectrum.sample_rate / 2.0)
  if total == 0.0 {
    return 0.0
  }
  let target = if fraction < 0.0 {
    0.0
  } else if fraction > 1.0 {
    1.0
  } else {
    fraction
  }
  let mut accumulated = 0.0
  for index = 0; index < spectrum.powers.length(); index = index + 1 {
    accumulated += spectrum.powers[index]
    if accumulated >= target * total {
      return spectrum.frequencies[index]
    }
  }
  spectrum.frequencies[spectrum.frequencies.length() - 1]
}

///|
/// Estimate the dominant frequency directly from a real signal.
pub fn dominant_frequency_of(
  data : Array[Double],
  sample_rate : Double,
) -> Double {
  dominant_frequency(hann_periodogram(data, sample_rate))
}

///|
/// Return the amplitude spectrum of a signal.
pub fn amplitude_spectrum(data : Array[Double]) -> Array[Double] {
  let bins = real_dft(data)
  let result = Array::new(capacity=bins.length())
  let scale = if data.is_empty() {
    0.0
  } else {
    2.0 / data.length().to_double()
  }
  for bin in bins {
    result.push(bin.magnitude() * scale)
  }
  result
}

///|
/// Compute spectral flux between adjacent frames.
pub fn spectral_flux(
  previous : Array[Double],
  current : Array[Double],
) -> Double {
  let old = amplitude_spectrum(previous)
  let new = amplitude_spectrum(current)
  let count = if old.length() < new.length() {
    old.length()
  } else {
    new.length()
  }
  let mut total = 0.0
  for index = 0; index < count; index = index + 1 {
    let difference = new[index] - old[index]
    if difference > 0.0 {
      total += difference * difference
    }
  }
  total.sqrt()
}

///|
/// Normalize a magnitude spectrum to unit maximum for feature comparison.
pub fn normalize_spectrum(amplitudes : Array[Double]) -> Array[Double] {
  let result = amplitudes.copy()
  let peak = maximum_sample(result)
  if peak == 0.0 {
    return result
  }
  for index = 0; index < result.length(); index = index + 1 {
    result[index] = result[index] / peak
  }
  result
}
