#! /bin/bash

cat requirements.txt requirements-dev.txt |while read p
do
  conda install -y "$p"
done

