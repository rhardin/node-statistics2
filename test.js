var stats = require('./lib/statistics2'),
    z = stats.pnormaldist(1 - (1 - 0.975) / 2),
    normal = stats.normaldist(0.27);

console.log(z);
console.log(normal); //0.60641987319804
