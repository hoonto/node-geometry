var g = require('./geometry.js');

var runtest = Module.cwrap('runtest', '', []);

runtest();

//var sphere_distance1 = Module.cwrap('sphere_distance1', 'number', ['number','number','number','number','number'])
//var toradians = 0.0174532925;
//var result = sphere_distance1(39.7391500 * toradians, -104.9847000 * toradians, 34.0522300 * toradians, -118.2436800 * toradians, 6371);
//console.log(result);

