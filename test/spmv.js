'use strict';

var chai = require('chai');
var assert = chai.assert;
// var blas1 = require('ndarray-blas-level1');
var RandMatGen = require('./util/rand-matrix-gen.js');
var ndarray = require('ndarray');
var assertCloseTo = require('./util/close-to');
var constants = require('./util/constants');

var gemv = require('../gemv');
var spmv = require('../spmv');

describe('SPMV (symmetric packed matrix vector product)', function () {
  var n = 4;
  var seed;
  var matGen = new RandMatGen(seed, Float64Array);
  var x = ndarray(new Float64Array(n), [n]);
  var y1 = ndarray(new Float64Array(n), [n]);
  var y2 = ndarray(new Float64Array(n), [n]);
  var len = ((n + 1) * n) / 2;
  var A = ndarray(new Float64Array(len), [len]);
  var B = ndarray(new Float64Array(n * n), [n, n]);

  it('upper-triangular SPMV', function () {
    for (var t = 0; t < constants.NUM_TESTS; ++t) {
      seed = matGen.setRandomSeed(36);
      matGen.makePackedMatrix(n, A);
      matGen.makeGeneralMatrix(1, n, x);
      matGen.makeSymmetricMatrixFromPacked(n, A, B);

      assert(spmv(A, n, x, y1, false));
      assert(gemv(1, B, x, 0, y2));
      assertCloseTo(y1, y2, constants.TEST_TOLERANCE, 'Failure seed value: "' + seed + '".');
    }
  });

  it('lower-triangular SPMV', function () {
    for (var t = 0; t < constants.NUM_TESTS; ++t) {
      seed = matGen.setRandomSeed(36);
      matGen.makePackedMatrix(n, A);
      matGen.makeGeneralMatrix(1, n, x);
      matGen.makeSymmetricMatrixFromPacked(n, A, B);

      assert(spmv(A, n, x, y1));
      assert(gemv(1, B, x, 0, y2));
      assertCloseTo(y1, y2, constants.TEST_TOLERANCE, 'Failure seed value: "' + seed + '".');
    }
  });
});
