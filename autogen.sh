#!/bin/sh

configure_args=""

while test $# -ge 1; do
    case $1 in
	--verbose) verbose=yes ;;
	--enable-*) configure_args="$configure_args $1" ;;
	--disable-*) configure_args="$configure_args $1" ;;
	--with-*) configure_args="$configure_args $1" ;;
	--without-*) configure_args="$configure_args $1" ;;
	*) echo "unknown argument $1"; exit 1 ;;
    esac
    shift
done

# paranoia: sometimes autoconf doesn't get things right the first time
autoreconf --verbose --install --symlink --force
autoreconf --verbose --install --symlink --force
autoreconf --verbose --install --symlink --force

config=good # hackery so darcs_test still outputs config.log w/failed configure

./configure --enable-maintainer-mode $configure_args || config=bad

if test x$verbose = xyes; then
    cat config.log
fi

test $config = bad && exit 1
