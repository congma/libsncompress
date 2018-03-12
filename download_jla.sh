#!/bin/sh
# Show usage.
usage ()
{
cat << _get_test_data_usage_EOF1317813562
usage: ${0} [-c] [-d directory] [-h]
  -h            show this help and exit
  -c            clean up after download
  -d directory  directory in which to write data files
                (default: testdata)
_get_test_data_usage_EOF1317813562
}


# Hard-coded SHA1 checksum of the complete list of required files.
# "Files with those SHA1 values present" is used as the invariant.
hardcoded_sha1 ()
{
cat << _get_test_data_checksums_EOF548586433
6e2df790acd5c23a27ce0fa53e5632cc12350d5b *covmat/C_bias.fits
521943edf2681f14e9ea246ed68c75e25860c8f2 *covmat/C_cal.fits
edfa73bf907f13a57ca61a4585ec02755354d5be *covmat/C_dust.fits
2fac326e3020bbecf51bd0d768d2e4c2c9080acc *covmat/C_host.fits
345f7434738dcf75ca7d2a8e4b0375703331c448 *covmat/C_model.fits
b161e42525e727807172774d8791d45d7d6a80e0 *covmat/C_nonia.fits
af9f74bcf6a07767845d492cce048525ae6cd538 *covmat/C_pecvel.fits
ab050d38e619fd1c83bbd5bb03a4d8dfdc1ac493 *covmat/C_stat.fits
6a7752cde98de993c278e2fda6015921eac1699b *jla_lcparams.txt
9db8919852c22dbf051b00b44883cc4392887947 *jla_mub.txt
21274ced39046ea5e897b265edefd816f475d6e4 *jla_mub_covmatrix.dat
_get_test_data_checksums_EOF548586433
}


# Newline-separated list of output file paths (relative to base dir).
output_files_in_basedir ()
{
    # NOTE: "cut" field starts from 1, not 0.
    hardcoded_sha1 | cut -d '*' -f 2-
}


# The invariant-check command in base dir.
invariant_in_basedir ()
{
    hardcoded_sha1 | "$SHA1SUMCMD" -c -
}


# Clean downloaded tarballs, if asked to do so.
clean_in_basedir ()
{
    if [ "$CLEAN" -ge 1 ]; then
	rm -rf covmat_v6.tgz jla_likelihood_v6.tgz
    fi
}


# Check if the tarball contains paths that may blow up when extracted.
has_sneaky_member ()
{
    tar ztf "$1" | grep -qE '^/|\.\.'
}


# Parse command-line arguments.
BASEDIR="testdata"
CLEAN=0
OPTIND=0

while getopts hcd: opt; do
    case $opt in
        h)
            usage
            exit 0
            ;;
        c)  CLEAN="$((CLEAN + 1))"
            ;;
        d)  BASEDIR="$OPTARG"
            ;;
        *)
            usage >&2
            exit 1
            ;;
    esac
done
shift "$((OPTIND - 1))"   # Discard the options and sentinel --


# Check if SHA1 calculator is useable.
if [ -z "$SHA1SUMCMD" ]; then
    if [ "$(uname)" = "Darwin" ]; then
	if which shasum > /dev/null 2>&1; then
	    # The shasum perl utility shipped with macos correctly fail on
	    # non-existent files, so it overrides.
	    SHA1SUMCMD="shasum"
	    # This is to avoid the bad implementation from homebrew's
	    # md5sha1sum installation.
	else
	    >&2 echo "Error: no useable SHA1 checksum utility present"
	    exit 1
	fi
    elif which sha1sum > /dev/null 2>&1; then
	SHA1SUMCMD="sha1sum"
    else
	>&2 echo "Error: no useable SHA1 checksum utility present"
	exit 1
    fi
fi


# Ensure that the directory exists
mkdir -p "$BASEDIR" > /dev/null 2>&1
if ! cd "$BASEDIR"; then
    >&2 echo "Error: cannot cd into the output directory"
    exit 1
fi


# Invariant already satisfied.
if { invariant_in_basedir > /dev/null 2>&1 ; }; then
    >&2 echo "Warning: files already present; will not download data tarballs"
    clean_in_basedir
    exit 0
fi


# Which command shall we use for downloading?
if [ -z "$DLCMD" ]; then
    if which wget > /dev/null 2>&1; then
	DLCMD="wget"
    elif which curl > /dev/null 2>&1; then
	DLCMD="curl"
    else
	>&2 echo "Error: no wget or curl present"
	exit 1
    fi
fi

__D=0
"$DLCMD" http://supernovae.in2p3.fr/sdss_snls_jla/covmat_v6.tgz
__D="$((__D + $?))"
sleep 2
"$DLCMD" http://supernovae.in2p3.fr/sdss_snls_jla/jla_likelihood_v6.tgz
__D="$((__D + $?))"

if [ "${__D}" != "0" ]; then
    >&2 echo "Error: download failed"
    exit "${__D}"
fi


# Check if downloaded files are tarbombs.
for bn in covmat_v6.tgz jla_likelihood_v6.tgz; do
    if has_sneaky_member "$bn"; then
	>&2 echo "Error: $bn appears sneaky"
	exit 1
    fi
done


# Extract tarballs
# Extract covmat tarball, overwriting existing files.
# (We already checked that the invariant hasn't been satisfied so far, so this
# is OK.)
# Get tarball members by filtering output file list already in the script.
output_files_in_basedir | grep 'covmat/C_.*\.fits' | \
    xargs tar -zxf covmat_v6.tgz
# Extract the text data files.  Also by filtering (and replacing).
output_files_in_basedir | grep 'jla_.*' | \
    sed -e 's#^#jla_likelihood_v6/data/#g' | \
    xargs tar --strip-components=2 -zxf jla_likelihood_v6.tgz
# Make extracted files read-only.
output_files_in_basedir | xargs chmod -w


# Check if the invariant is satisfied.
if ! { invariant_in_basedir > /dev/null ; }; then
    >&2 echo "Error: some downloaded and extracted files do not match checksum"
    __ERR=1
else
    __ERR=0
fi


# Cleanup and exit
clean_in_basedir
exit "${__ERR}"
