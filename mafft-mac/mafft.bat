#! /bin/sh

pushd "`dirname "$0"`" > /dev/null 2>&1; rootdir="$PWD"; export rootdir; popd > /dev/null 2>&1;
MAFFT_BINARIES="$rootdir"/mafftdir/libexec; export MAFFT_BINARIES;

arguments=""
while [ $# -gt 1 ];
do
	arguments=$arguments" "$1
	shift
done

"$rootdir"/mafftdir/bin/mafft $arguments "$1"
# $1 can have space in file name
