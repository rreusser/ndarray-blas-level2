'use strict';

var blas1 = require('ndarray-blas-level1');

module.exports = sbmv;

function sbmv (A, k, x, y, fromLower, alpha, beta) {
  var n = A.shape[0];

  if (k > n - 1) {
    // number of superdiagonals cannot exceed matrix dimensions.
    return false;
  }
  if (n !== A.shape[1]) {
    // matrix must be symmetric
    return false;
  }
  if (n !== x.shape[0]) {
    // x dimension must agree with A
    return false;
  }
  if (n !== y.shape[0]) {
    // y dimension must agree with A
    return false;
  }
  var lower = fromLower || true;
  var alpha0 = alpha === undefined ? 1 : alpha;
  var beta0 = beta === undefined ? 0 : beta;

  var i = 0;
  var j = 0;
  if (beta0 === 0) {
    for (i = 0; i < y.shape[0]; ++i) {
      y.set(i, 0);
    }
  } else if (beta0 !== 1) {
    blas1.scal(beta0, y);
  }
  if (alpha0 === 0) {
    return true;
  }

  var temp1 = 0.0;
  var temp2 = 0.0;
  var l = 0.0;
  var kplus1 = k + 1;
  if (lower) {
    for (j = 0; j < n; ++j) {
      temp1 = alpha0 * x.get(j);
      temp2 = 0;
      y.set(j, y.get(j) + temp1 * A.get(0, j));
      l = 1 - j;
      for (i = j + 1; i < Math.min(n, j + k); ++i) {
        y.set(i, y.get(i) + temp1 * A.get(l + i, j));
        temp2 = temp2 + A.get(l + i, j) * x.get(i);
      }
      y.set(j, y.get(j) + alpha0 * temp2);
    }
  } else {
    for (j = 0; j < n; ++j) {
      temp1 = alpha0 * x.get(j);
      temp2 = 0;
      l = kplus1 - j;
      for (i = Math.max(l, j - k); i < j - 1; ++i) {
        y.set(i, y.get(i) + temp1 * A.get(l + i, j));
        temp2 = temp2 + A.get(l + i, j) * x.get(i);
      }
      y.set(j, y.get(j) + temp1 * A.get(kplus1, j) + alpha0 * temp2);
    }
  }

  return true;
}
