# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.1] - 2023-06-13

### Changed
- Move to pyproject.toml packaging for compatibility with newer PIP versions
- Fixed incorrect parameter name in gamtools call_windows

## [2.0.0] - 2022
### Added
- Allow the doit dependencies database to be reset
- Calculate coverage on individual bamfiles instead of all at once
- Calculate intrinsic bias of a contact matrix
- Added functions for calculating normalized pointwise mutual information (NPMI)
- Added functions for checking the quality of a segregation table
- Added a CONTRIBUTORS file

### Changed
- Changed default matrix type from Dprime to NPMI
- Use base-pair instead of read-count coverage for window calling
- Use Orphan windows algorithm to call positive windows

## [1.1.1] - 2019-02-08
### Added
- Added a CHANGELOG

### Changed
- Fix problem processing uncompressed fastas
- Fix bug where pandas converts sample names to numeric types

### Removed
- Removed command line option for MACS window calling
- Removed support for Python3.4

## [1.1.0] - 2017-05-06
### Added
- Added fixed threshold option to process_nps
- Added functions for extracting sub-regions of a matrix

### Changed
- Enable coveralls
- Enable codeclimate


## [1.0.0] - 2017-03-04
### Added
- First release of GAMtools package
- Created a master command line script
- Added support for gzipped input/output
- Added a command to calculate enrichments
- Added a command to convert different matrix formats
- Added a command for circular permutation
- Added a command for getting compaction
- Added a command for getting radial position
- Added QC scripts and automatic exclusion of bad samples
- Added documentation
- Added unit tests
- Added tutorial
- Added continuous integration through tox and Travis

### Changed
- Made code PEP8 compliant
- Made code python3 compatible

### Removed
- Old command line scripts

[1.0.0]: https://github.com/pombo-lab/gamtools/compare/nature-2017...v1.0.0
[1.1.0]: https://github.com/pombo-lab/gamtools/compare/v1.0.0...v1.1.0
[1.1.1]: https://github.com/pombo-lab/gamtools/compare/v1.1.0...v1.1.1

