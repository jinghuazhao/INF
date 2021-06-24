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

mkdocs build
mkdocs gh-deploy

git add .gitignore
git commit -m ".gitignore"
git add README.md
git commit -m "README"
git add index.html
git commit -m "index.html"
git add docs
git commit -m "docs"
git add mkdocs.yml
git commit -m "mkdocs.yml"
git push
