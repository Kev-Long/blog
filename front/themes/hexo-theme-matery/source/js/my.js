
const yaml = require('node-yaml')

console.log("here")

yaml.read('../../languagues/default.yml', (err, data) => {
    const res = JSON.parse(JSON.stringify(data));
    res.postCategoryTitle = '';
    yaml.write('../../languagues/default.yml', res, (err) => {
        if(err) {
        console.log(err.toString());
        }
    });
    });