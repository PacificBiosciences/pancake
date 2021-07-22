V:=1.4.21
R:=racon-v$V
B:=$R/build-meson
#export LDFLAGS?="-static-libstdc++ -static-libgcc"
#F:=--default-library static

all:
	${MAKE} untar
	rm -rf $B
	${MAKE} conf
	${MAKE} install
	rm -rf $B
	rm -rf $R
	rm -rf $R.tar.gz
untar: $R
$R: $R.tar.gz
	tar xvfz $<
$R.tar.gz:
	curl -OL https://github.com/lbcb-sci/racon/releases/download/$V/racon-v$V.tar.gz
conf:
	mkdir -p $B && cd $B && meson setup --buildtype=release -Dc_args=-O3 -Dtests=true --prefix=${PREFIX_ARG} . ..
install:
	ninja -v -C $B install
