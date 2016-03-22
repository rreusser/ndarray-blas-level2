'use strict';

module.exports.gemv = require('./gemv.js');
module.exports.gbmv = require('./gbmv.js');
module.exports.symv = require('./symv.js');
module.exports.sbmv = require('./sbmv.js');
module.exports.spmv = require('./spmv.js');
module.exports.trmv = require('./trmv.js');
module.exports.tbmv = require('./tbmv.js');
module.exports.trsv = require('./trsv.js');
module.exports.tbsv = require('./tbsv.js');
module.exports.tpsv = require('./tpsv.js');
module.exports.ger = require('./ger.js');
module.exports.syr = require('./syr.js');
module.exports.spr = require('./spr.js');
module.exports.syr2 = require('./syr2.js');
module.exports.spr2 = require('./spr2.js');
module.exports.trmv_lower = function (A, x) {
  console.warn('trmv_lower is deprecated. Please use the \'isLower\' flag with trmv.');
  return module.exports.trmv(A, x, true);
};
module.exports.trsv_lower = function (A, x) {
  console.warn('trsv_lower is deprecated. Please use the \'isLower\' flag with trsv.');
  return module.exports.trsv(A, x, true);
};
