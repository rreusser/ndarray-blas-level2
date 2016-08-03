'use strict';

var blas1 = require('ndarray-blas-level1');

module.exports = spmv;

function spmv (A, n, x, y, fromLower, alpha, beta) {
  var lower = fromLower || true;
  var alpha0 = alpha === undefined ? 1 : alpha;
  var beta0 = beta === undefined ? 0 : beta;

  if (A.shape[0] !== n * (n + 1) / 2) {
    throw new Error('Packed matrix does not match arguments.');
  }
  if (n !== x.shape[0]) {
    throw new Error('x dimension must agree with A.');
  }
  if (n !== y.shape[0]) {
    throw new Error('y dimension must agree with A.');
  }
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

  var kk = 0;
  var temp1 = 0;
  var temp2 = 0;
  var k = 0;
  if (lower) {
    for (j = 0; j < n; ++j) {
      temp1 = alpha0 * x.get(j);
      temp2 = 0;
      y.set(j, y.get(j) + temp1 * A.get(kk));
      k = kk + 1;
      for (i = j + 1; i < n; ++i) {
        y.set(i, y.get(i) + temp1 * A.get(k));
        temp2 = temp2 + A.get(k) * x.get(i);
        k = k + 1;
      }
      y.set(j, y.get(j) + alpha0 * temp2);
      kk = kk + (n - j);
    }
  } else {
    for (j = 0; j < n; ++j) {
      temp1 = alpha0 * x.get(j);
      temp2 = 0;
      k = kk;
      for (i = 0; i < j - 1; ++i) {
        y.set(i, y.get(i) + temp1 * A.get(k));
        temp2 = temp2 + A.get(k) * x.get(i);
        k = k + 1;
      }
      y.set(j, y.get(j) + (alpha0 * temp2) + (temp1 * A.get(kk + j)));
      kk = kk + j + 1;
    }
  }
  return true;
}
