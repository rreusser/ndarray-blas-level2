'use strict';

var RandMatGen = require('./util/rand-matrix-gen.js');
var ndarray = require('ndarray');
var assertCloseTo = require('./util/close-to');
var constants = require('./util/constants');

var gemv = require('../gemv');
var sbmv = require('../sbmv');

describe('SBMV (symmetric banded matrix-vector product)', function () {
  var n = 15;
  var k = 3;
  var seed;
  var matGen = new RandMatGen(seed, Float64Array);
  var x = ndarray(new Float64Array(n), [n]);
  var y1 = ndarray(new Float64Array(n), [n]);
  var y2 = ndarray(new Float64Array(n), [n]);
  var A = ndarray(new Float64Array(n * n), [n, n]);

  it('upper-triangular sbmv', function () {
    for (var t = 0; t < constants.NUM_TESTS; ++t) {
      seed = matGen.setRandomSeed(36);
      matGen.makeSymmBandedMatrix(n, k, A);

      sbmv(A, k, x, y1, false);
      gemv(1, A, x, 0, y2);
      assertCloseTo(y1, y2, constants.TEST_TOLERANCE, 'Failure seed value: "' + seed + '".');
    }
  });

  it('lower-triangular sbmv', function () {
    for (var t = 0; t < constants.NUM_TESTS; ++t) {
      seed = matGen.setRandomSeed(36);
      matGen.makeSymmBandedMatrix(n, k, A);

      sbmv(A, k, x, y1);
      gemv(1, A, x, 0, y2);
      assertCloseTo(y1, y2, constants.TEST_TOLERANCE, 'Failure seed value: "' + seed + '".');
    }
  });
});
