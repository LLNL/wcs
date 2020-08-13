#!/bin/bash

if [ $# -ne 1 ] ; then
  working_dir=.
else
  working_dir=$1
fi

pushd ${working_dir} > /dev/null

if [ ! -d iheap ] || [ ! -f iheap/iheap.h ] ; then
  git clone https://github.com/yangle/iheap.git
  patch --no-backup-if-mismatch iheap/iheap.h < patch_iheap_h.txt
fi

popd > /dev/null
