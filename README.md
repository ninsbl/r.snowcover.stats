# r.snowcover.stats: GRASS GIS module example in Python

[![GRASS GIS module](https://img.shields.io/badge/GRASS%20GIS-module-%23009000)](https://grass.osgeo.org/)
[![CI](https://github.com/ninsbl/r.snowcover.stats/workflows/CI/badge.svg)](https://github.com/ninsbl/r.snowcover.stats/actions?query=workflow%3A%22CI%22)
[![Code quality check](https://github.com/ninsbl/r.snowcover.stats/workflows/Code%20quality%20check/badge.svg)](https://github.com/ninsbl/r.snowcover.stats/actions?query=workflow%3A%22Code%20quality%20check%22)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Black code style check](https://github.com/ninsbl/r.snowcover.stats/workflows/Black%20code%20style%20check/badge.svg)](https://github.com/ninsbl/r.snowcover.stats/actions?query=workflow%3A%22Black%20code%20style%20check%22)
[![Deploy online documentation](https://github.com/ninsbl/r.snowcover.stats/workflows/Deploy%20online%20documentation/badge.svg)](https://github.com/ninsbl/r.snowcover.stats/actions?query=workflow%3A%22Deploy%20online%20documentation%22)

This is an example of a Python script to run in actinia using GRASS GIS parser and deployment

Documentation of the script is generated in a github action and available here:
https://ninsbl.github.io/r.snowcover.stats/

### Getting the GitHub Actions work

If you have the repository on GitHub, you can also reuse the GitHub
Actions defined in the repository (under `.github`). Initially, most of
them will fail, but once you do the renaming, most of them should start
working.

For the workflow uploading documentation to GitHub Pages to
work, you will need you to
[set up Deploy key and a Secret](https://github.com/marketplace/actions/github-pages-action#1-add-ssh-deploy-key)
for your repository. Once the keys are in place, the online documentation
will be published as GitHub Pages website automatically.
The URL for the website is available in the Settings of your repository.
