///|
/// Boundary coverage for statistics, windows, and numerical helpers.
test "advanced statistics boundaries" {
  assert_true(mean_samples([]) == 0.0, "empty mean")
  assert_true(variance_samples([]) == 0.0, "empty variance")
  assert_true(
    quantile_samples([1.0, 2.0, 3.0], -1.0) == 1.0,
    "low quantile clamp",
  )
  assert_true(
    quantile_samples([1.0, 2.0, 3.0], 2.0) == 3.0,
    "high quantile clamp",
  )
  assert_near(median_samples([3.0, 1.0, 2.0]), 2.0, 0.0000001, "median")
  assert_true(zero_crossing_count([0.0]) == 0, "single crossing")
  assert_true(
    normalize_peak([0.0, 0.0]) == [0.0, 0.0],
    "constant normalization",
  )
  assert_true(
    scale_samples([2.0, 2.0], 4.0, 8.0) == [4.0, 4.0],
    "constant scaling",
  )
  let stats = describe_samples([-2.0, 0.0, 2.0])
  assert_true(stats.count == 3, "stats count")
  assert_near(stats.peak_to_peak, 4.0, 0.0000001, "stats range")
}

///|
test "window and resampling boundaries" {
  assert_true(window_count(0, 4, 2, false) == 0, "empty window count")
  assert_true(window_count(5, 4, 2, false) == 1, "full window count")
  assert_true(window_count(5, 4, 2, true) == 3, "partial window count")
  assert_true(
    extract_window([1.0, 2.0], -1, 4, 0.0) == [0.0, 1.0, 2.0, 0.0],
    "padded window",
  )
  assert_true(hann_window(0).length() == 0, "empty hann")
  assert_true(hann_window(1) == [1.0], "single hann")
  assert_near(
    cumulative_integral([1.0, 1.0, 1.0], 1.0)[2],
    2.0,
    0.0000001,
    "integral",
  )
  assert_near(
    differentiate_samples([0.0, 2.0, 4.0], 1.0)[1],
    2.0,
    0.0000001,
    "derivative",
  )
  assert_true(resample_linear([], 10.0, 20.0).is_empty(), "empty resample")
  assert_true(
    decimate_samples([1.0, 2.0, 3.0], 0).is_empty(),
    "invalid decimation",
  )
}

///|
test "spectral and filter response" {
  let wave = sine_wave(64, 4.0, 32.0, 1.0, 0.0)
  let spectrum = hann_periodogram(wave, 32.0)
  assert_true(dominant_frequency(spectrum) > 2.0, "dominant frequency positive")
  assert_true(spectral_entropy(spectrum) >= 0.0, "entropy non-negative")
  assert_true(spectral_rolloff(spectrum, 0.0) >= 0.0, "rolloff lower bound")
  let lowpass_result = design_fir_lowpass(32.0, 5.0, 16, fir_hann())
  match lowpass_result {
    Ok(filter) => {
      assert_true(filter.coefficients.length() % 2 == 1, "odd FIR taps")
      assert_true(
        apply_fir(wave, filter).length() == wave.length(),
        "FIR length",
      )
    }
    Err(message) => abort(message)
  }
  match design_fir_lowpass(32.0, 0.0, 7, fir_hann()) {
    Ok(_) => abort("invalid FIR accepted")
    Err(_) => ()
  }
  let noisy = [0.0, 0.0, 10.0, 0.0, 0.0]
  assert_true(median_filter(noisy, 3)[2] == 0.0, "median spike removal")
  assert_true(gaussian_kernel(2, 1.0).length() == 5, "gaussian kernel size")
}

///|
test "quality and validation boundaries" {
  let empty = Trace::new("", "", 1.0, 0.0, [])
  let report = assess_quality(empty, 1.0, 0.0)
  assert_true(!is_good_quality(report), "empty quality finding")
  let invalid_rate = Trace::new("X", "Y", 0.0, 0.0, [1.0])
  match validate_trace(invalid_rate) {
    Ok(_) => abort("invalid rate accepted")
    Err(_) => ()
  }
  let findings = validate_samples([0.0, 2.0], 1.0, 1.0)
  assert_true(validation_warning_count(findings) == 1, "range warning")
  assert_true(validate_sac_roundtrip_size(empty), "empty SAC size invariant")
  assert_true(
    validate_miniseed_roundtrip_size(Trace::new("A", "B", 10.0, 0.0, [1.0])),
    "miniSEED size",
  )
}

///|
test "event detection and synthetic boundaries" {
  let event = SyntheticEvent::new(20, 3.0, 3.0, 20, 0.0)
  let trace = synthetic_trace("STA", "BHZ", 64, 32.0, event, 0.001)
  let config = default_event_config(trace.delta(), trace.start_time)
  let events = detect_events(trace.data, config)
  assert_true(events.length() >= 0, "event detector completes")
  assert_true(
    sum_samples(impulse_signal(8, 100, 1.0, 1)) == 0.0,
    "out of range impulse",
  )
  assert_true(step_signal(4, 10, 1.0, 0.0) == [0.0, 0.0, 0.0, 0.0], "late step")
  assert_true(first_arrival_index([], 1.0, 2) is None, "empty arrival")
  let candidate = strongest_event(events)
  match candidate {
    Some(value) => assert_true(event_duration(value) >= 0.0, "event duration")
    None => ()
  }
}

///|
test "travel time and model boundaries" {
  let model = iaspei91()
  assert_true(validate_velocity_model(model) is Ok(_), "IASPEI model valid")
  assert_true(
    validate_velocity_model(
      VelocityModel::new("bad", [Layer::new(2.0, 3.0, 4.0)]),
    )
    is Err(_),
    "invalid velocity model",
  )
  assert_near(
    travel_time_steps(model, 0.0, 0.0, wave_phase_p(), 1),
    0.0,
    0.0000001,
    "zero path",
  )
  assert_true(
    vertical_travel_time(model, 20.0, wave_phase_p()) > 0.0,
    "vertical time",
  )
  assert_true(
    ray_path_points(10.0, 20.0, 1).length() == 2,
    "minimum ray points",
  )
  let curve = sample_travel_time_curve(
    model,
    10.0,
    0.0,
    20.0,
    4,
    wave_phase_p(),
  )
  assert_true(curve.samples.length() == 4, "travel curve samples")
}

///|
test "format and binary metadata boundaries" {
  match parse_ascii_samples("# header\n1.0  2.0\n3.0\n") {
    Ok(values) => assert_true(values == [1.0, 2.0, 3.0], "ASCII parsing")
    Err(message) => abort(message)
  }
  match parse_ascii_samples("bad\n") {
    Ok(_) => abort("bad ASCII accepted")
    Err(_) => ()
  }
  let trace = Trace::new("ST", "BHZ", 20.0, 2.0, [0.1, -0.2, 0.3])
  let bytes = to_sac_bytes(trace, true)
  match parse_sac_header(bytes) {
    Ok(header) => {
      assert_true(header.sample_count == 3, "SAC header count")
      assert_true(header.endian is SacBigEndian, "SAC endian")
    }
    Err(message) => abort(message)
  }
  let mini = to_miniseed_bytes(trace, 7, true)
  match parse_miniseed_header(mini) {
    Ok(header) => assert_true(header.sample_count == 3, "miniSEED header count")
    Err(message) => abort(message)
  }
  match
    parse_time_amplitude_csv("time,amplitude\n0,1\n0.5,2\n1,3\n", "A", "B") {
    Ok(parsed) => assert_near(parsed.sampling_rate, 2.0, 0.0000001, "CSV rate")
    Err(message) => abort(message)
  }
}

///|
test "correlation, streaming, and catalog" {
  let accumulator = RollingAccumulator::new()
  accumulator.push(1.0)
  accumulator.push(3.0)
  assert_near(accumulator.mean(), 2.0, 0.0000001, "online mean")
  assert_near(accumulator.variance(), 1.0, 0.0000001, "online variance")
  let detector = StreamingStaLta::new(2, 4)
  ignore(detector.push(1.0))
  ignore(detector.push(2.0))
  assert_true(detector.push(4.0) > 0.0, "streaming ratio")
  let ring = SampleRing::new(2)
  ring.push(1.0)
  ring.push(2.0)
  ring.push(3.0)
  assert_true(ring.snapshot() == [2.0, 3.0], "ring eviction")
  let correlation = cross_correlate([0.0, 1.0, 0.0], [1.0], 2)
  assert_true(correlation.coefficients.length() == 5, "correlation bins")
  let catalog = EventCatalog::new()
  ignore(catalog.add(1.0, 0.0, 0.0, 5.0, 2.0, 3, 0.8))
  assert_true(catalog.above_magnitude(1.5).length() == 1, "catalog filter")
  assert_true(catalog.mean_magnitude() == 2.0, "catalog mean")
}

///|
test "phase, response, and pipeline boundaries" {
  let trace = Trace::new(
    "P",
    "BHZ",
    20.0,
    0.0,
    sine_wave(64, 2.0, 20.0, 1.0, 0.0),
  )
  assert_true(
    instantaneous_phase(trace.data).length() == trace.npts(),
    "phase length",
  )
  assert_true(instantaneous_frequency([], 20.0).is_empty(), "empty frequency")
  assert_true(estimate_period(trace.data, 20.0) >= 0.0, "period estimate")
  let response = generic_broadband_response()
  assert_true(validate_response(response) is Ok(_), "response valid")
  assert_true(
    calibrate_trace(trace, response, 1.0).npts() == trace.npts(),
    "calibration length",
  )
  match analyze_trace(trace, default_processing_config()) {
    Ok(result) =>
      assert_true(result.processed.npts() == trace.npts(), "pipeline length")
    Err(message) => abort(message)
  }
  let workload = BenchmarkWorkload::new("small", 32, 16.0, 1, 0.01)
  let summary = run_benchmark_workload(workload, 1)
  assert_true(summary.samples_processed == 32, "benchmark samples")
}

///|
test "geographic and comparison boundaries" {
  let a = StationPosition::new("A", 0.0, 0.0, 0.0)
  let b = StationPosition::new("B", 0.0, 1.0, 0.0)
  assert_true(surface_distance_km(a, a) == 0.0, "zero geographic distance")
  assert_true(surface_distance_km(a, b) > 100.0, "geographic distance")
  assert_true(azimuth_degrees(a, b) >= 0.0, "azimuth range")
  let left = Trace::new("A", "X", 10.0, 0.0, [1.0, 2.0, 3.0])
  let right = Trace::new("A", "X", 10.0, 0.0, [1.0, 2.0, 3.0])
  match compare_traces(left, right, 1) {
    Ok(comparison) =>
      assert_true(comparison.rms_error == 0.0, "identical trace comparison")
    Err(message) => abort(message)
  }
  assert_true(agreement_ratio([1.0], [1.0], 0.0) == 1.0, "agreement")
}
