# agc-istream
Header-only C++ library for std::istream compatible access to
[agc](https://github.com/refresh-bio/agc) compressed genome archives
(work in progress).

## Dependencies
- [agc](https://github.com/refresh-bio/agc) >= v3.0
- [zstd](https://github.com/facebook/zstd)

## Usage
Opean the agc archive `archive.agc` and read the sample `sample` into an agc::istream with
```
agc::istream stream("archive.agc");
stream.find("sample");
```
After initialization, the contents of `agc::istream` can be accessed in
any manner that is supported by the `std::istringstream` class. If the
sample is found in the archive, agc::istream status is set to fail.

## License
The source code from this project is subject to the terms of the MIT
license. A copy of the MIT license is supplied with the project, or
can be obtained at https://opensource.org/licenses/MIT.
