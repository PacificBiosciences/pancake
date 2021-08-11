# pancake

A new PacBio HiFi overlapper. Work in progress.

## Dependencies

```
Gcc / Clang with C++17 support
Meson >=0.52.0
Ninja
Boost
```

Note: to use the GPU functionalities, please specify the `-Dgpu-cuda=true` flag when configuring the Meson build.
```
ENABLED_GPU_CUDA=true bash -vex scripts/ci/configure_with_fallback.sh
make build
```
Alternatively, the following Make rule can be used:
```
make conf-gpu
make build
```

To successfully compile the CUDA-related code (GenomeWorks), a modified version of Meson is needed. It can be obtained and set up like this:
```
git clone https://github.com/SoapGentoo/meson.git -b cuda-fixes
cd meson
ln -sf meson.py meson
export PATH=$(pwd):$PATH
```

## Instructions

```
export CC=/usr/bin/gcc-7
export CXX=/usr/bin/g++-7
make all
```

## FAQ
This project uses boost.ICL which has some known incompatibility with gcc c++17 support. If you run into the following error -

```
/usr/include/boost/icl/type_traits/type_to_string.hpp:56:12: error: partial specialization of ‘struct boost::icl::type_to_string<Unary<Type> >’ after instantiation of ‘struct boost::icl::type_to_string<std::__cxx11::basic_string<char> >’ [-fpermissive]
   56 |     struct type_to_string<Unary<Type> >
      |            ^~~~~~~~~~~~~~~~~~~~~~~~~~~~
```

track down which project build the error is coming from (e.g. top level project, or any of the subprojects) and update the `meson.build` file to include `'-fno-new-ttp-matching'` in the `CXXFLAGS` section.
