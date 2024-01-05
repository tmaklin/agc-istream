# agcstr
(work in progress) Read samples from [agc](https://github.com/refresh-bio/agc) archives using C++ iostreams.

## Usage
Opean the agc archive `archive.agc` and read the sample `sample` into an agc_istream buffer with
```
agc_istream stream("archive.agc");
stream.find("sample");
```
After initialization, the contents of `agc_istream` can be accessed in
any manner that is supported by the `std::istringstream` class.
