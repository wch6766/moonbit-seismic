///| Constant PI for waveform generation in tests.
const PI : Double = 3.141592653589793

///| Warning-free assertion helper for boolean conditions.
fn assert_true(cond : Bool, msg : String) -> Unit {
  if !cond {
    abort("Assertion failed: " + msg)
  }
}

///| Warning-free assertion helper for comparing doubles within epsilon.
fn assert_near(val1 : Double, val2 : Double, eps : Double, msg : String) -> Unit {
  let diff = (val1 - val2).abs()
  if diff > eps {
    abort("Assertion failed: " + msg + " (expected " + val2.to_string() + ", got " + val1.to_string() + ", diff " + diff.to_string() + ")")
  }
}

test "Trace basic stats and preprocessing" {
  let data = [1.0, 2.0, 3.0, 4.0, 5.0]
  let trace = Trace::new("STA01", "BHZ", 10.0, 0.0, data)
  
  assert_true(trace.npts() == 5, "npts")
  assert_near(trace.delta(), 0.1, 0.0000001, "delta")
  assert_near(trace.duration(), 0.4, 0.0000001, "duration")
  
  assert_near(trace.mean(), 3.0, 0.0000001, "mean")
  assert_near(trace.std(), 1.58113883, 0.000001, "std")
  
  // Test demean
  trace.demean()
  assert_near(trace.mean(), 0.0, 0.0000001, "demean mean")
  assert_near(trace.data[0], -2.0, 0.0000001, "demean sample 0")
  assert_near(trace.data[4], 2.0, 0.0000001, "demean sample 4")
  
  // Test detrend
  let trend_data = [5.0, 7.0, 9.0, 11.0, 13.0]
  let trend_trace = Trace::new("STA01", "BHZ", 1.0, 0.0, trend_data)
  trend_trace.detrend()
  for i = 0; i < 5; i = i + 1 {
    assert_near(trend_trace.data[i], 0.0, 0.0000001, "detrend sample " + i.to_string())
  }
}

test "Interpolation algorithms" {
  let x = [0.0, 1.0, 2.0]
  let y = [10.0, 20.0, 40.0]
  
  // Linear interpolation
  let r1 = linear_interpolate(x, y, 0.5)
  match r1 {
    Ok(val) => assert_near(val, 15.0, 0.0000001, "linear interpolate 0.5")
    Err(e) => abort(e)
  }
  let r2 = linear_interpolate(x, y, 1.5)
  match r2 {
    Ok(val) => assert_near(val, 30.0, 0.0000001, "linear interpolate 1.5")
    Err(e) => abort(e)
  }
  let r_left = linear_interpolate(x, y, -1.0)
  match r_left {
    Ok(val) => assert_near(val, 10.0, 0.0000001, "linear interpolate left extrapolation")
    Err(e) => abort(e)
  }
  
  // Cubic Spline
  let spline_res = CubicSpline::new(x, y)
  match spline_res {
    Ok(spline) => {
      // Natural boundary conditions mean c[0] and c[2] should be 0.0
      assert_near(spline.c[0], 0.0, 0.0000001, "spline natural boundary c[0]")
      assert_near(spline.c[2], 0.0, 0.0000001, "spline natural boundary c[2]")
      
      // Interpolated values
      let val_05 = spline.evaluate(0.5)
      let val_15 = spline.evaluate(1.5)
      
      assert_near(spline.evaluate(0.0), 10.0, 0.0000001, "spline evaluate 0.0")
      assert_near(spline.evaluate(1.0), 20.0, 0.0000001, "spline evaluate 1.0")
      assert_near(spline.evaluate(2.0), 40.0, 0.0000001, "spline evaluate 2.0")
      
      assert_true(val_05 > 10.0 && val_05 < 20.0, "spline evaluate 0.5 bounds")
      assert_true(val_15 > 20.0 && val_15 < 40.0, "spline evaluate 1.5 bounds")
    }
    Err(e) => abort(e)
  }
}

test "Travel-time curves" {
  let model = iaspei91()
  assert_true(model.name == "IASPEI91", "model name")
  
  let tp = travel_time(model, 10.0, 100.0, false)
  let ts = travel_time(model, 10.0, 100.0, true)
  
  assert_true(tp > 0.0, "P travel time positive")
  assert_true(ts > tp, "S travel time greater than P")
  assert_near(tp, 100.5 / 5.86, 2.0, "P travel time estimate")
}

test "Waveform filtering" {
  let n = 200
  let sampling_rate = 50.0
  let data = Array::new(capacity=n)
  
  for i = 0; i < n; i = i + 1 {
    let t = i.to_double() / sampling_rate
    let val = @math.sin(2.0 * PI * 1.0 * t) + 0.5 * @math.sin(2.0 * PI * 15.0 * t)
    data.push(val)
  }
  
  let trace = Trace::new("FILT", "BHZ", sampling_rate, 0.0, data)
  lowpass(trace, 3.0, true)
  
  assert_true(trace.npts() == 200, "filtered trace length")
}

test "STA/LTA picking and triggering" {
  let n = 200
  let data = Array::make(n, 0.05)
  for i = 100; i < n; i = i + 1 {
    data[i] = 1.0 * @math.sin((i - 100).to_double() * 0.2)
  }
  
  let ratio = classic_sta_lta(data, 5, 20)
  
  assert_true(ratio[50] < 2.0, "pre-event ratio")
  assert_true(ratio[110] > 2.0, "event ratio")
  
  let triggers = trigger_onset(ratio, 2.0, 1.2, 0.01, 0.0)
  assert_true(triggers.length() >= 1, "trigger count")
  assert_true(triggers[0].on_index >= 95 && triggers[0].on_index <= 110, "trigger onset index")
}

test "Cross-correlation lag detection" {
  let sig = [0.0, 0.0, 1.0, 2.0, 1.0, 0.0, 0.0]
  let temp = [1.0, 2.0, 1.0]
  
  let res = cross_correlate(sig, temp, 4)
  assert_true(res.max_lag == 2, "max lag")
  assert_near(res.max_coefficient, 1.0, 0.0000001, "max coefficient")
}

test "SAC file parsing and writing roundtrip" {
  let data = [0.1, 0.2, -0.3, 0.4, 0.5]
  let orig_trace = Trace::new("STN01", "BHE", 20.0, 10.0, data)
  
  // 1. Little endian roundtrip
  let sac_bytes_le = to_sac_bytes(orig_trace, false)
  let parsed_res_le = parse_sac(sac_bytes_le)
  match parsed_res_le {
    Ok(trace) => {
      assert_true(trace.station == "STN01", "station le")
      assert_true(trace.channel == "BHE", "channel le")
      assert_near(trace.sampling_rate, 20.0, 0.00001, "sampling rate le")
      assert_near(trace.start_time, 10.0, 0.00001, "start time le")
      assert_true(trace.npts() == 5, "npts le")
      assert_near(trace.data[2], -0.3, 0.00001, "data sample 2 le")
    }
    Err(e) => abort(e)
  }
  
  // 2. Big endian roundtrip
  let sac_bytes_be = to_sac_bytes(orig_trace, true)
  let parsed_res_be = parse_sac(sac_bytes_be)
  match parsed_res_be {
    Ok(trace) => {
      assert_true(trace.station == "STN01", "station be")
      assert_true(trace.channel == "BHE", "channel be")
      assert_near(trace.sampling_rate, 20.0, 0.00001, "sampling rate be")
      assert_near(trace.start_time, 10.0, 0.00001, "start time be")
      assert_true(trace.npts() == 5, "npts be")
      assert_near(trace.data[2], -0.3, 0.00001, "data sample 2 be")
    }
    Err(e) => abort(e)
  }
}

test "miniSEED file parsing and writing roundtrip" {
  let data = [10.0, 20.0, 30.0, 40.0, 50.0]
  let orig_trace = Trace::new("ST02", "BHZ", 100.0, 0.0, data)
  
  let mseed_bytes = to_miniseed_bytes(orig_trace, 1, false)
  let parsed_res = parse_miniseed_record(mseed_bytes, false)
  match parsed_res {
    Ok(trace) => {
      assert_true(trace.station == "ST02", "mseed station")
      assert_true(trace.channel == "BHZ", "mseed channel")
      assert_near(trace.sampling_rate, 100.0, 0.00001, "mseed sampling rate")
      assert_true(trace.npts() == 5, "mseed npts")
      assert_near(trace.data[3], 40.0, 0.00001, "mseed data sample 3")
    }
    Err(e) => abort(e)
  }
}
