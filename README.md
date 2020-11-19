# pancake

A new PacBio HiFi overlapper. Work in progress.

## Dependencies

```
meson
ninja
```

## Instructions

```
export CC=/usr/bin/gcc-7
export CXX=/usr/bin/g++-7
make build
```

## FAQ
This project uses boost.ICL which has some known incompatibility with gcc c++17 support. If you run into the following error - 

```
/usr/include/boost/icl/type_traits/type_to_string.hpp:56:12: error: partial specialization of ‘struct boost::icl::type_to_string<Unary<Type> >’ after instantiation of ‘struct boost::icl::type_to_string<std::__cxx11::basic_string<char> >’ [-fpermissive]
   56 |     struct type_to_string<Unary<Type> >
      |            ^~~~~~~~~~~~~~~~~~~~~~~~~~~~
```

track down which project build the error is coming from (e.g. top level project, or any of the subprojects) and update the `meson.build` file to include `'-fno-new-ttp-matching'` in the `CXXFLAGS` section.
