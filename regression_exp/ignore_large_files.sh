#!/bin/bash

# Find all files > 100MB and add to .gitignore
find . -type f -size +100M -not -path "./.git/*" | sed 's|^\./||' >> .gitignore

# Sort and remove duplicates
sort -u -o .gitignore .gitignore

echo "Added all files >100MB to .gitignore. You can now safely commit changes!"
