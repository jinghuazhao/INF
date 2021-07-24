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
pip install mkdocs-mermaid2-plugin
}

# setup

mkdocs build
mkdocs gh-deploy

git add .gitignore
git commit -m ".gitignore"
git add README.md
git commit -m "README"
git add docs
git commit -m "docs"
git add mkdocs.yml
git commit -m "mkdocs.yml"
git push
