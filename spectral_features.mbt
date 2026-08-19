///|
/// Features commonly used to classify or triage seismic windows.
pub struct SpectralFeatureSet {
  dominant_frequency : Double
  centroid : Double
  spread : Double
  entropy : Double
  rolloff_85 : Double
  low_band_power : Double
  mid_band_power : Double
  high_band_power : Double
  spectral_flatness : Double
}

///|
/// Return a frequency bin index for a requested frequency.
pub fn frequency_bin(spectrum : Spectrum, frequency : Double) -> Int {
  if spectrum.frequencies.is_empty() || frequency <= 0.0 {
    return 0
  }
  let estimate = (frequency *
  spectrum.input_length.to_double() /
  spectrum.sample_rate).to_int()
  if estimate < 0 {
    0
  } else if estimate >= spectrum.frequencies.length() {
    spectrum.frequencies.length() - 1
  } else {
    estimate
  }
}

///|
/// Return a geometric mean divided by arithmetic mean spectral flatness.
pub fn spectral_flatness(spectrum : Spectrum) -> Double {
  if spectrum.powers.is_empty() {
    return 0.0
  }
  let mut logarithmic = 0.0
  let mut arithmetic = 0.0
  let mut count = 0
  for power in spectrum.powers {
    let safe = if power <= 0.0 { 1.0e-30 } else { power }
    logarithmic += @math.ln(safe)
    arithmetic += power
    count += 1
  }
  if arithmetic <= 0.0 || count == 0 {
    0.0
  } else {
    @math.exp(logarithmic / count.to_double()) /
    (arithmetic / count.to_double())
  }
}

///|
/// Extract a stable feature vector from a one-sided spectrum.
pub fn extract_spectral_features(spectrum : Spectrum) -> SpectralFeatureSet {
  let nyquist = spectrum.sample_rate / 2.0
  {
    dominant_frequency: dominant_frequency(spectrum),
    centroid: spectral_centroid(spectrum),
    spread: spectral_spread(spectrum),
    entropy: spectral_entropy(spectrum),
    rolloff_85: spectral_rolloff(spectrum, 0.85),
    low_band_power: band_power(spectrum, 0.0, nyquist / 3.0),
    mid_band_power: band_power(spectrum, nyquist / 3.0, 2.0 * nyquist / 3.0),
    high_band_power: band_power(spectrum, 2.0 * nyquist / 3.0, nyquist),
    spectral_flatness: spectral_flatness(spectrum),
  }
}

///|
/// Extract spectral features directly from a signal.
pub fn features_from_samples(
  data : Array[Double],
  sample_rate : Double,
) -> SpectralFeatureSet {
  extract_spectral_features(hann_periodogram(data, sample_rate))
}

///|
/// Return local maxima in an amplitude spectrum above a threshold.
pub fn spectral_peaks(
  spectrum : Spectrum,
  minimum_amplitude : Double,
) -> Array[(Double, Double)] {
  let result = Array::new()
  if spectrum.amplitudes.length() < 3 {
    return result
  }
  for index = 1; index < spectrum.amplitudes.length() - 1; index = index + 1 {
    let amplitude = spectrum.amplitudes[index]
    if amplitude >= minimum_amplitude &&
      amplitude >= spectrum.amplitudes[index - 1] &&
      amplitude >= spectrum.amplitudes[index + 1] {
      result.push((spectrum.frequencies[index], amplitude))
    }
  }
  result
}

///|
/// Estimate harmonicity as the ratio of the strongest peak to total power.
pub fn harmonicity(spectrum : Spectrum) -> Double {
  let peaks = spectral_peaks(spectrum, 0.0)
  let mut strongest = 0.0
  for pair in peaks {
    let (_, amplitude) = pair
    if amplitude > strongest {
      strongest = amplitude
    }
  }
  let total = band_power(spectrum, 0.0, spectrum.sample_rate / 2.0).sqrt()
  if total == 0.0 {
    0.0
  } else {
    strongest / total
  }
}

///|
/// Interpolate a magnitude at a frequency between adjacent bins.
pub fn interpolate_spectral_amplitude(
  spectrum : Spectrum,
  frequency : Double,
) -> Double {
  if spectrum.frequencies.is_empty() {
    return 0.0
  }
  if frequency <= spectrum.frequencies[0] {
    return spectrum.amplitudes[0]
  }
  let last = spectrum.frequencies.length() - 1
  if frequency >= spectrum.frequencies[last] {
    return spectrum.amplitudes[last]
  }
  let index = frequency_bin(spectrum, frequency)
  let left = if spectrum.frequencies[index] > frequency {
    index - 1
  } else {
    index
  }
  let right = if left + 1 >= spectrum.frequencies.length() {
    left
  } else {
    left + 1
  }
  let width = spectrum.frequencies[right] - spectrum.frequencies[left]
  if width == 0.0 {
    spectrum.amplitudes[left]
  } else {
    let fraction = (frequency - spectrum.frequencies[left]) / width
    spectrum.amplitudes[left] +
    fraction * (spectrum.amplitudes[right] - spectrum.amplitudes[left])
  }
}

///|
/// Return the fraction of power inside a closed band.
pub fn band_power_fraction(
  spectrum : Spectrum,
  low_frequency : Double,
  high_frequency : Double,
) -> Double {
  let total = band_power(spectrum, 0.0, spectrum.sample_rate / 2.0)
  if total == 0.0 {
    0.0
  } else {
    band_power(spectrum, low_frequency, high_frequency) / total
  }
}

///|
/// Return true when the dominant bin is inside the requested band.
pub fn dominant_in_band(
  spectrum : Spectrum,
  low_frequency : Double,
  high_frequency : Double,
) -> Bool {
  let frequency = dominant_frequency(spectrum)
  frequency >= low_frequency && frequency <= high_frequency
}
