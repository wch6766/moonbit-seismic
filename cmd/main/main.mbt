///| Entry point for the CLI or runnable program in this template.
fn main {
  println("==================================================")
  println("   MoonBit Seismic Wave Signal Processing Demo    ")
  println("==================================================")
  
  // 1. Create a synthetic seismic trace with a transient wave
  let n = 200
  let sampling_rate = 50.0
  let data = Array::new(capacity=n)
  
  // Background noise + event starting at sample 100
  for i = 0; i < n; i = i + 1 {
    let t = i.to_double() / sampling_rate
    let noise = 0.05 * @math.sin(2.0 * 3.141592653589793 * 15.0 * t) // high frequency noise
    let event = if i >= 100 {
      1.0 * @math.sin((i - 100).to_double() * 0.1)
    } else {
      0.0
    }
    data.push(noise + event)
  }
  
  let trace = @seismic.Trace::new("ST01", "BHZ", sampling_rate, 0.0, data)
  println("1. Created synthetic Trace:")
  println("   Station: " + trace.station + ", Channel: " + trace.channel)
  println("   Samples: " + trace.npts().to_string() + ", Sampling Rate: " + trace.sampling_rate.to_string() + " Hz")
  println("   Mean: " + trace.mean().to_string() + ", Std Dev: " + trace.std().to_string())
  
  // 2. Filter the trace to remove noise
  let filtered = trace.clone()
  @seismic.lowpass(filtered, 5.0, true)
  println("\n2. Applied 2nd-order Lowpass Filter (Cutoff: 5.0 Hz, Zero-phase):")
  println("   Original Std Dev: " + trace.std().to_string())
  println("   Filtered Std Dev: " + filtered.std().to_string())
  
  // 3. STA/LTA pick arrival
  let ratio = @seismic.classic_sta_lta(filtered.data, 5, 20)
  let triggers = @seismic.trigger_onset(ratio, 2.0, 1.2, filtered.delta(), filtered.start_time)
  println("\n3. Ran STA/LTA Arrival Picker:")
  if triggers.length() > 0 {
    // Access fields from Trigger struct via getters or directly
    // since fields are public by default in structural matches or package visibility
    // Let's print trigger details using a loop
    for t in triggers {
      // Print first trigger info
      println("   Detected Event Arrival!")
      println("   Onset Time: " + t.on_time.to_string() + " s")
      println("   Offset Time: " + t.off_time.to_string() + " s")
      break
    }
  } else {
    println("   No trigger found.")
  }
  
  // 4. Magnitude Estimation
  // Suppose peak amplitude measured is 5.5 mm on Wood-Anderson seismograph at 120 km epicentral distance
  let ml = @seismic.local_magnitude(5.5, 120.0)
  println("\n4. Earthquake Magnitude Estimation:")
  match ml {
    Ok(val) => println("   Richter Local Magnitude (M_L): " + val.to_string())
    Err(err) => println("   Error: " + err)
  }
  
  // 5. Travel-time curves for IASPEI91 velocity model
  let model = @seismic.iaspei91()
  let tp = @seismic.travel_time(model, 15.0, 80.0, false)
  let ts = @seismic.travel_time(model, 15.0, 80.0, true)
  println("\n5. 1D Travel-Time Calculations (IASPEI91 model, source depth 15km, distance 80km):")
  println("   P-wave travel time: " + tp.to_string() + " s")
  println("   S-wave travel time: " + ts.to_string() + " s")
  
  // 6. Velocity Model Interpolation
  let depths = [0.0, 10.0, 20.0, 35.0]
  let velocities = [5.5, 6.0, 6.5, 8.0]
  let spline_res = @seismic.CubicSpline::new(depths, velocities)
  println("\n6. Velocity Model Interpolation:")
  match spline_res {
    Ok(spline) => {
      let v_at_15 = spline.evaluate(15.0)
      println("   Cubic Spline interpolated velocity at 15.0 km depth: " + v_at_15.to_string() + " km/s")
    }
    Err(e) => println("   Spline Error: " + e)
  }
  
  // 7. SAC Serialization
  let sac_bytes = @seismic.to_sac_bytes(trace, false)
  println("\n7. SAC File Serialization:")
  println("   Successfully serialized trace to SAC binary format (" + sac_bytes.length().to_string() + " bytes).")
  
  println("==================================================")
}
