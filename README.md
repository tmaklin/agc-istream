# agcstr
Header-only C++ library for std::istream compatible access to
[agc](https://github.com/refresh-bio/agc) compressed genome archives
(work in progress) .

## Usage
Opean the agc archive `archive.agc` and read the sample `sample` into an agc::istream with
```
agc::istream stream("archive.agc");
stream.find("sample");
```
After initialization, the contents of `agc_istream` can be accessed in
any manner that is supported by the `std::istringstream` class. If the
sample is found in the archive, agc::istream status is set to fail.
