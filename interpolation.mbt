///| Perform 1D linear interpolation.
/// Extrapolates to boundary values if x_target is out of bounds.
pub fn linear_interpolate(
  x : Array[Double],
  y : Array[Double],
  x_target : Double
) -> Result[Double, String] {
  let n = x.length()
  if n != y.length() {
    return Err("Input arrays x and y must have the same length")
  }
  if n == 0 {
    return Err("Arrays must not be empty")
  }
  
  if x_target <= x[0] {
    return Ok(y[0])
  }
  if x_target >= x[n - 1] {
    return Ok(y[n - 1])
  }
  
  // Binary search to find the containing interval
  let mut low = 0
  let mut high = n - 1
  while low < high - 1 {
    let mid = (low + high) / 2
    if x[mid] <= x_target {
      low = mid
    } else {
      high = mid
    }
  }
  
  let x0 = x[low]
  let x1 = x[high]
  let y0 = y[low]
  let y1 = y[high]
  let t = (x_target - x0) / (x1 - x0)
  Ok(y0 + t * (y1 - y0))
}

///| Cubic Spline Interpolator struct.
pub struct CubicSpline {
  x : Array[Double]
  a : Array[Double]
  b : Array[Double]
  c : Array[Double]
  d : Array[Double]
}

///| Construct a Natural Cubic Spline from data points x and y.
/// Points in x must be strictly increasing.
pub fn CubicSpline::new(x : Array[Double], y : Array[Double]) -> Result[CubicSpline, String] {
  let n = x.length()
  if n != y.length() {
    return Err("x and y must have the same length")
  }
  if n < 3 {
    return Err("Cubic spline requires at least 3 points")
  }
  
  // Ensure x is strictly increasing
  for i = 1; i < n; i = i + 1 {
    if x[i] <= x[i - 1] {
      return Err("x array must be strictly increasing")
    }
  }
  
  let h = Array::make(n - 1, 0.0)
  for i = 0; i < n - 1; i = i + 1 {
    h[i] = x[i + 1] - x[i]
  }
  
  // Vector of equation targets
  let alpha = Array::make(n - 1, 0.0)
  for i = 1; i < n - 1; i = i + 1 {
    alpha[i] = 3.0 / h[i] * (y[i + 1] - y[i]) - 3.0 / h[i - 1] * (y[i] - y[i - 1])
  }
  
  // Solve tridiagonal system using Thomas Algorithm
  let l = Array::make(n, 1.0)
  let mu = Array::make(n, 0.0)
  let z = Array::make(n, 0.0)
  
  for i = 1; i < n - 1; i = i + 1 {
    l[i] = 2.0 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1]
    mu[i] = h[i] / l[i]
    z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i]
  }
  
  let c = Array::make(n, 0.0)
  let b = Array::make(n - 1, 0.0)
  let d = Array::make(n - 1, 0.0)
  
  // Back substitution
  let mut j = n - 2
  while j >= 0 {
    c[j] = z[j] - mu[j] * c[j + 1]
    b[j] = (y[j + 1] - y[j]) / h[j] - h[j] * (c[j + 1] + 2.0 * c[j]) / 3.0
    d[j] = (c[j + 1] - c[j]) / (3.0 * h[j])
    if j == 0 {
      break
    }
    j = j - 1
  }
  
  let a = Array::new(capacity=n)
  for v in y {
    a.push(v)
  }
  
  Ok({ x, a, b, c, d })
}

///| Evaluate the Cubic Spline at a specific target.
pub fn CubicSpline::evaluate(self : CubicSpline, x_target : Double) -> Double {
  let n = self.x.length()
  if x_target <= self.x[0] {
    return self.a[0]
  }
  if x_target >= self.x[n - 1] {
    return self.a[n - 1]
  }
  
  // Binary search to find containing interval
  let mut low = 0
  let mut high = n - 1
  while low < high - 1 {
    let mid = (low + high) / 2
    if self.x[mid] <= x_target {
      low = mid
    } else {
      high = mid
    }
  }
  
  let dx = x_target - self.x[low]
  self.a[low] + self.b[low] * dx + self.c[low] * dx * dx + self.d[low] * dx * dx * dx
}
