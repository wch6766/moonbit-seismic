///| Represents a detected seismic event window.
pub struct Trigger {
  on_index : Int
  off_index : Int
  on_time : Double
  off_time : Double
}

///| Compute the classic STA/LTA (Short-Term Average / Long-Term Average) ratio
/// using an efficient O(N) sliding-window recursive algorithm.
pub fn classic_sta_lta(
  data : Array[Double],
  n_sta : Int,
  n_lta : Int
) -> Array[Double] {
  let n = data.length()
  let ratio = Array::make(n, 1.0)
  
  if n == 0 || n_sta <= 0 || n_lta <= 0 || n_lta <= n_sta {
    return ratio
  }
  
  let mut sum_sta = 0.0
  let mut sum_lta = 0.0
  
  // Initialize sums for the first n_lta samples
  for i = 0; i < n_lta; i = i + 1 {
    let val = data[i].abs()
    sum_lta = sum_lta + val
    if i >= n_lta - n_sta {
      sum_sta = sum_sta + val
    }
  }
  
  let sta_denom = n_sta.to_double()
  let lta_denom = n_lta.to_double()
  
  let initial_lta = sum_lta / lta_denom
  ratio[n_lta - 1] = if initial_lta == 0.0 { 1.0 } else { (sum_sta / sta_denom) / initial_lta }
  
  // Slide the window across the remaining samples
  for i = n_lta; i < n; i = i + 1 {
    let val_in = data[i].abs()
    
    // Slide LTA window
    let val_out_lta = data[i - n_lta].abs()
    sum_lta = sum_lta + val_in - val_out_lta
    
    // Slide STA window
    let val_out_sta = data[i - n_sta].abs()
    sum_sta = sum_sta + val_in - val_out_sta
    
    let lta = sum_lta / lta_denom
    let sta = sum_sta / sta_denom
    
    ratio[i] = if lta == 0.0 { 1.0 } else { sta / lta }
  }
  
  ratio
}

///| Detect event trigger windows based on onset/offset thresholds.
pub fn trigger_onset(
  ratio : Array[Double],
  on_threshold : Double,
  off_threshold : Double,
  dt : Double,
  start_time : Double
) -> Array[Trigger] {
  let n = ratio.length()
  let triggers = Array::new()
  let mut active = false
  let mut on_idx = 0
  
  for i = 0; i < n; i = i + 1 {
    let r = ratio[i]
    if !active && r >= on_threshold {
      active = true
      on_idx = i
    } else if active && r < off_threshold {
      active = false
      let off_idx = i
      let on_time = start_time + on_idx.to_double() * dt
      let off_time = start_time + off_idx.to_double() * dt
      triggers.push({
        on_index: on_idx,
        off_index: off_idx,
        on_time,
        off_time,
      })
    }
  }
  
  // Close the final trigger if it is still active at the end of the trace
  if active {
    let off_idx = n - 1
    let on_time = start_time + on_idx.to_double() * dt
    let off_time = start_time + off_idx.to_double() * dt
    triggers.push({
      on_index: on_idx,
      off_index: off_idx,
      on_time,
      off_time,
    })
  }
  
  triggers
}
