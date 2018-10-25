#!/bin/bash

set -e -o pipefail

mkdir -p $INSTALL_PATH/ngenebio_ngs_etc
if [ ! -f $INSTALL_PATH/ngenebio_ngs_etc/ngenebio_ngs ]; then
	ln -s $NGENEBIO_NGS_PATH $INSTALL_PATH/ngenebio_ngs_etc/ngenebio_ngs
fi

sync

if [[ "$1" == "minimal" ]]; then
	conda install --name env --override-channels -y \
		-q -c bioconda -c conda-forge -c defaults -c r \
		--file "$NGENEBIO_NGS_PATH/requirements.txt"
else
	conda install --name env --override-channels -y \
		-q -c bioconda -c conda-forge -c defaults -c r \
		--file "$NGENEBIO_NGS_PATH/requirements-conda-tests.txt"
fi

conda clean -y --all
