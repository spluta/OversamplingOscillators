// #include <cmath>

// int kaiser( float beta, int M, float *window )
// {
//   // # Docstring adapted from NumPy's kaiser function
//   // if _len_guards(M):
//   //     return np.ones(M)
//   // M, needs_trunc = _extend(M, sym)

//   int result = true;

//   // n = np.arange(0, M)
//   float n[M];
//   for ( int i = 0; i < M; i++ )
//     n[i] = i;
//   // alpha = (M - 1) / 2.0
//   float alpha = (M - 1) / 2.0;
//   // w = (special.i0(beta * np.sqrt(1 - ((n - alpha) / alpha) ** 2.0)) /
//   //      special.i0(beta))
//   for ( int i = 0; i < M; i++ ) {
//     float p = pow( (n[i] - alpha) / alpha, 2 );
//     window[i] = i0( beta * sqrt(1 - p) ) / i0( beta );
//   }

//   return result;
// }