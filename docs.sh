#!/usr/bin/bash

function nonweb ()
{
git add doc
git commit -m "Documents and auxiliary files"
git add cardio
git commit -m "Work on CEU's cardio"
git add csd3
git commit -m "Work on CSD3"
git add rsid
git commit -m "Final work on CSD3 with rsids in the downstream analysis"
git add tryggve
git commit -m "Notes and programs on TRYGGVE"
git push
}

function setup()
{
module load python/3.7
source ~/COVID-19/py37/bin/activate
}

function install()
{
  pip install mkdocs-mermaid2-plugin
  pip install mkdocs-section-index
}

setup

mkdocs build
mkdocs gh-deploy

for f in .gitignore README.md docs mkdocs.yml
do
  git add ${f}
  git commit -m "${f}"
done
git push
