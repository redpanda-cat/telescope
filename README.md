# Telescope

Telescope is a dashboard Jinja template generator to show components from [Planetarium](https://github.com/shahcompbio/planetarium).

## How to Use

You can specify what dashboard to use via `App.js`

### Developer

```
yarn install
```

```
yarn start
```

### Building

First compile all the code - this should end with files in `/build`

```
yarn build
```

Then we will want to merge all assets into one file - this will be compiled with the given name in `dist/<name>`

This is done with a python script, so ensure you have the right environment with `requirements.txt` added.

```
python build_template.py <name of dashboard> 
```


## Current Dashboards

- [Sankey](https://github.com/shahcompbio/sankey_clone)
